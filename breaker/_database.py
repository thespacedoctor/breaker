#!/usr/local/bin/python
# encoding: utf-8
"""
*Class to setup database object for the breaker*

:Author:
    David Young

:Date Created:
    October 29, 2015
"""
################# GLOBAL IMPORTS ####################
import sys
import os
os.environ['TERM'] = 'vt100'
import readline
import glob
import pickle
import time
import MySQLdb as ms
from docopt import docopt
from fundamentals import tools, times


class database():

    """
    *The worker class for the database module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
    """
    # INITIALISATION

    def __init__(
            self,
            log,
            settings=False,

    ):
        self.log = log
        log.debug("instansiating a new '_database' object")
        self.settings = settings
        return None

    def close(self):
        del self
        return None

    # METHOD ATTRIBUTES
    def get(self):
        """
        *get the database object*

        **Return:**
            - ``self.transientsDbConn, self.ps1gwDbConn, self.cataloguesDbConn`` -- three database connections
        """
        self.log.debug('starting the ``get`` method')
        self._setup_database_connections()
        self.log.debug('completed the ``get`` method')
        return self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn

    def _setup_database_connections(
            self):
        """
        *setup database connections for transient and catalogue databases*
        """
        self.log.debug('starting the ``_setup_database_connections`` method')

        from subprocess import Popen, PIPE, STDOUT

        tunnelDatabases = {}
        for db in self.settings["database settings"]:
            port = self.settings["database settings"][db]["port"]
            if "tunnel" in str(port):
                if port not in tunnelDatabases.keys():
                    tunnelDatabases[port] = []
                tunnelDatabases[port].append(
                    db)

        # SETUP TUNNEL IF REQUIRED
        if "ssh tunnels" in self.settings:
            tunnels = self.settings["ssh tunnels"]
            for tunnelName in tunnels:
                tunnel = self.settings["ssh tunnels"][tunnelName]
                # TEST TUNNEL DOES NOT ALREADY EXIST
                sshPort = tunnel["port"]
                for db in tunnelDatabases[tunnelName]:
                    connected = self._checkServer(
                        self.settings["database settings"][db]["host"], sshPort)
                    if connected:
                        break
                if connected:
                    self.log.debug('ssh tunnel already exists - moving on')
                else:
                    # GRAB TUNNEL SETTINGS FROM SETTINGS FILE
                    ru = tunnel["remote user"]
                    rip = tunnel["remote ip"]
                    rh = tunnel["remote datbase host"]

                    cmd = "ssh -fnN %(ru)s@%(rip)s -L %(sshPort)s:%(rh)s:3306" % locals()
                    p = Popen(cmd, shell=True, close_fds=True)
                    output = p.communicate()[0]
                    self.log.debug('output: %(output)s' % locals())

                    # TEST CONNECTION - QUIT AFTER SO MANY TRIES
                    connected = False
                    count = 0
                    while not connected:
                        for db in tunnelDatabases[tunnelName]:
                            connected = self._checkServer(
                                self.settings["database settings"][db]["host"], sshPort)
                            if connected:
                                break
                        time.sleep(1)
                        count += 1
                        if count == 5:
                            self.log.error(
                                'cound not setup tunnel to remote datbase' % locals())
                            sys.exit(0)

        # SETUP A DATABASE CONNECTION FOR THE ps1gw
        host = self.settings["database settings"][
            "ps1gw"]["host"]
        user = self.settings["database settings"][
            "ps1gw"]["user"]
        passwd = self.settings["database settings"][
            "ps1gw"]["password"]
        dbName = self.settings["database settings"][
            "ps1gw"]["db"]
        port = self.settings["database settings"][
            "ps1gw"]["port"]
        if "tunnel" in str(port):
            port = self.settings["ssh tunnels"][port]["port"]
        thisConn = ms.connect(
            host=host,
            user=user,
            passwd=passwd,
            db=dbName,
            port=port,
            use_unicode=True,
            charset='utf8'
        )
        thisConn.autocommit(True)
        self.log.debug('ps1gwDbConn: %s' % (thisConn,))
        self.ps1gwDbConn = thisConn

        # SETUP DATABASE CONNECTION FOR WAVE DATABASE
        host = self.settings["database settings"][
            "ligo_virgo_waves"]["host"]
        user = self.settings["database settings"][
            "ligo_virgo_waves"]["user"]
        passwd = self.settings["database settings"][
            "ligo_virgo_waves"]["password"]
        dbName = self.settings["database settings"][
            "ligo_virgo_waves"]["db"]
        port = self.settings["database settings"][
            "ligo_virgo_waves"]["port"]
        if "tunnel" in str(port):
            port = self.settings["ssh tunnels"][port]["port"]
        thisConn = ms.connect(
            host=host,
            user=user,
            passwd=passwd,
            db=dbName,
            port=port,
            use_unicode=True,
            charset='utf8'
        )
        thisConn.autocommit(True)
        self.log.debug('ligo_virgo_wavesDbConn: %s' % (thisConn,))
        self.ligo_virgo_wavesDbConn = thisConn

        # SETUP DATABASE CONNECTION FOR CATALOGUE DATABASE
        host = self.settings["database settings"][
            "catalogues"]["host"]
        user = self.settings["database settings"][
            "catalogues"]["user"]
        passwd = self.settings["database settings"][
            "catalogues"]["password"]
        dbName = self.settings["database settings"][
            "catalogues"]["db"]
        port = self.settings["database settings"][
            "catalogues"]["port"]
        if "tunnel" in str(port):
            port = self.settings["ssh tunnels"][port]["port"]
        thisConn = ms.connect(
            host=host,
            user=user,
            passwd=passwd,
            db=dbName,
            port=port,
            use_unicode=True,
            charset='utf8'
        )
        thisConn.autocommit(True)
        self.log.debug('catalogues database connection: %s' % (thisConn,))
        self.cataloguesDbConn = thisConn

        self.log.debug('completed the ``_setup_database_connections`` method')
        return None

    def _checkServer(self, address, port):
        """
        *Check that the TCP Port we've decided to use for tunnelling is available*
        """
        self.log.debug('starting the ``_checkServer`` method')

        # CREATE A TCP SOCKET
        import socket
        s = socket.socket()
        self.log.debug(
            """Attempting to connect to `%(address)s` on port `%(port)s`""" % locals())
        try:
            s.connect((address, port))
            self.log.debug(
                """Connected to `%(address)s` on port `%(port)s`""" % locals())
            return True
        except socket.error, e:
            self.log.warning(
                """Connection to `%(address)s` on port `%(port)s` failed - try again: %(e)s""" % locals())
            return False

        return None

    # xt-class-method


if __name__ == '__main__':
    main()
