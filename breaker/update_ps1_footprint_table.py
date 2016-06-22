#!/usr/local/bin/python
# encoding: utf-8
"""
*Update the PS1 footprint table in breaker database and associate with GWs*

:Author:
    David Young

:Date Created:
    October 29, 2015
"""
################# GLOBAL IMPORTS ####################
import sys
import os
os.environ['TERM'] = 'vt100'
import numpy as np
from fundamentals import tools, times
from fundamentals.mysql import readquery, insert_list_of_dictionaries_into_database_tables, writequery
from astrocalc.coords import unit_conversion
from HMpTy.mysql import add_htm_ids_to_mysql_database_table


class update_ps1_footprint_table():
    """
    *Update the PS1 footprint table in breaker database and associate with GWs*

    Metadata for each GW event should be found in the settings file and are used when associating the telescope pointings in the database with a GW event. For example, here are the metadata for the first GW burst:

    .. code-block:: yaml 

        gravitational waves:
            G184098:
                time:
                    mjdStart: 57279.90
                    mjdEnd: 57369.90
                plot:
                    raRange: 48.  # CENTRAL WIDTH IN DEGREES
                    decRange: 45.  # CENTRAL HEIGHT IN DEGREES
                    centralCoordinate: [141., 0.0]
                mapPath: "/Users/Dave/Dropbox/notes/astronotes-wiki/projects/GW-G184098/Alternate-Skymap-Stats-For-G184098/maps/LALInference_skymap.fits"

    Both sky-area and the time-range are used in the gw-pointing associations.

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``updateNed`` -- do you want to update NED stream in database (can take a LONG time). Default *False*

    **Usage:**

        .. code-block:: python 

            from breaker import update_ps1_footprint_table
            dbUpdater = update_ps1_footprint_table(
                log=log, 
                settings=settings,
                updateNed=False
            )
            dbUpdater.get()
    """

    def __init__(
            self,
            log,
            settings=False,
            updateNed=False
    ):
        self.log = log
        log.debug("instansiating a new 'update_ps1_footprint_table' object")
        self.settings = settings
        self.updateNed = updateNed

        # xt-self-arg-tmpx

        # SETUP THE VARIOUS DATABASE CONNECTIONS REQUIRED
        from breaker import database
        db = database(
            log=self.log,
            settings=self.settings
        )
        self.ligo_virgo_wavesDbConn, self.ps1gwDbConn, self.cataloguesDbConn = db.get()

        return None

    def get(self):
        """
        *Import the new PS1 pointings and query NED for new data*

        This method:

            * Imports the new PS1 pointings from the PS1 database, 
            * attempts to label these pointings with the ID for an associated GW, 
            * queries NED for new data covered by the sky-area of these pointings and 
            * adds data to the NED stream database table

        See the ``update_ps1_footprint_table`` class of usage info.

        **Return:**
            - None
        """
        self.log.info('starting the ``get`` method')

        self.import_new_ps1_pointings()
        self.label_pointings_with_gw_ids()
        self.populate_ps1_subdisk_table()
        if self.updateNed:
            self.update_ned_database_table()

        self.log.info('completed the ``get`` method')
        return None

    def import_new_ps1_pointings(
            self):
        """
        *Import any new PS1 GW pointings from the ps1gw database into the ``ps1_pointings`` table of the Ligo-Virgo Waves database*

        **Return:**
            - None

         **Usage:**

            .. code-block:: python 

                # IMPORT NEW PS1 POINTINGS FROM PS1 GW DATABASE INTO LIGO-VIRGO WAVES DATABASE
                from breaker import update_ps1_footprint_table
                dbUpdater = update_ps1_footprint_table(
                    log=log,
                    settings=settings
                )
                dbUpdater.import_new_ps1_pointings()
        """
        self.log.info('starting the ``import_new_ps1_pointings`` method')

        # SELECT ALL OF THE POINTING INFO REQUIRED FROM THE ps1gw DATABASE
        sqlQuery = u"""
            select imageid, m.exptime exp_time, truncate(mjd_obs, 8) mjd, m.fpa_ra as ra, m.fpa_dec as decl, fpa_filter
              from tcs_cmf_metadata m
              where m.fpa_ra != "NaN" and m.fpa_dec != "NaN"
            group by exp_time, mjd, fpa_object, fpa_comment, fpa_ra, fpa_dec, fpa_filter;
        """ % locals()
        rows = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ps1gwDbConn,
            quiet=False
        )

        # TIDY RESULTS BEFORE IMPORT
        entries = []

        converter = unit_conversion(
            log=self.log
        )
        for row in rows:
            e = {}
            e["raDeg"] = converter.ra_sexegesimal_to_decimal(
                ra=row["ra"]
            )
            e["decDeg"] = converter.dec_sexegesimal_to_decimal(
                dec=row["decl"]
            )
            e["exp_time"] = row["exp_time"]
            e["mjd"] = row["mjd"]
            e["filter"] = row["fpa_filter"][0]
            e["ps1_exp_id"] = row["imageid"]
            entries.append(e)

        # ADD THE NEW RESULTS TO THE ps1_pointings TABLE
        insert_list_of_dictionaries_into_database_tables(
            dbConn=self.ligo_virgo_wavesDbConn,
            log=self.log,
            dictList=entries,
            dbTableName="ps1_pointings",
            uniqueKeyList=["raDeg", "decDeg", "mjd"],
            dateModified=False,
            batchSize=2500
        )

        # APPEND HTMIDs TO THE ps1_pointings TABLE
        add_htm_ids_to_mysql_database_table(
            raColName="raDeg",
            declColName="decDeg",
            tableName="ps1_pointings",
            dbConn=self.ligo_virgo_wavesDbConn,
            log=self.log,
            primaryIdColumnName="primaryId"
        )

        print "PS1 pointings synced between `tcs_cmf_metadata` and `ps1_pointings` database tables"

        self.log.info('completed the ``import_new_ps1_pointings`` method')
        return None

    def label_pointings_with_gw_ids(
            self):
        """
        *Attempt to label the PS1 pointing with the GW IDs*

        The GW metadata used to associate PS1 pointings is taken from the settings file

        **Return:**
            - None

         **Usage:**

            .. code-block:: python 

                # ATTEMPT TO LABEL PS1 POINTINGS IN DATABASE WITH A GW ID 
                from breaker import update_ps1_footprint_table
                dbUpdater = update_ps1_footprint_table(
                    log=log,
                    settings=settings
                )
                dbUpdater.label_pointings_with_gw_ids()
        """
        self.log.info('starting the ``label_pointings_with_gw_ids`` method')

        # WAVE METADATA FOUND IN SETTINGS FILE
        for wave in self.settings["gravitational waves"]:

            # UNPACK THE PLOT PARAMETERS FROM THE SETTINGS FILE
            centralCoordinate = self.settings["gravitational waves"][
                wave]["plot"]["centralCoordinate"]
            raRange = self.settings["gravitational waves"][
                wave]["plot"]["raRange"]
            decRange = self.settings["gravitational waves"][
                wave]["plot"]["decRange"]

            raMax = (centralCoordinate[0] + raRange / 2.) + 5.
            raMin = (centralCoordinate[0] - raRange / 2.) - 5.
            decMax = (centralCoordinate[1] + decRange / 2.) + 5.
            decMin = (centralCoordinate[1] - decRange / 2.) - 5.

            mjdLower = self.settings["gravitational waves"][
                wave]["time"]["mjdStart"]
            mjdUpper = self.settings["gravitational waves"][
                wave]["time"]["mjdEnd"]

            if raMin > 0. and raMax < 360.:
                raWhere = """(raDeg > %(raMin)s and raDeg < %(raMax)s)""" % locals(
                )
            elif raMin < 0.:
                raMin2 = raMin + 360.
                raWhere = """((raDeg > 0. and raDeg < %(raMax)s) or raDeg > %(raMin2)s)""" % locals(
                )
            elif raMax > 360.:
                raMax2 = raMax - 360.
                raWhere = """((raDeg > %(raMin)s and raDeg < 360.) or raDeg < %(raMax2)s)""" % locals(
                )

            decWhere = """(decDeg > %(decMin)s and  decDeg < %(decMax)s)""" % locals(
            )

            mjdWhere = "(mjd>%(mjdLower)s and mjd<%(mjdUpper)s)" % locals()

            sqlQuery = u"""
                update ps1_pointings set gw_id = "%(wave)s" where %(raWhere)s and %(decWhere)s and %(mjdWhere)s and gw_id is null
            """ % locals()
            writequery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
            )

        sqlQuery = u"""
            select count(*) as count from ps1_pointings where gw_id is null;
        """ % locals()

        count = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn,
            quiet=False
        )[0]["count"]

        print "PS1 pointings labelled with their associated GW id"

        if count == 0:
            print "    Note all pointings have been labelled with GW ID"
        else:
            print "    %(count)s pointings remain unlabelled with a GW ID" % locals()

        self.log.info('completed the ``label_pointings_with_gw_ids`` method')
        return None

    def populate_ps1_subdisk_table(
            self):
        """
        *Calculate 49 subdisks for each of the PS1 pointings (used to query NED in manageable sized batches) and add them to the ``ps1_pointings_subdisks`` table of the database*

        .. image:: http://i.imgur.com/y3G0aax.png
            :width: 600 px

        **Return:**
            - None

         **Usage:**

            .. code-block:: python 

                # SPLIT PS1 POINTINGS INTO SUB-DISKS AND ADD TO LV DATABASE 
                from breaker import update_ps1_footprint_table
                dbUpdater = update_ps1_footprint_table(
                    log=log,
                    settings=settings
                )
                dbUpdater.populate_ps1_subdisk_table()
        """
        self.log.info(
            'starting the ``populate_ps1_subdisk_table`` method')

        # SELECT THE PS1 POINTINGS NEEDING SUBDISKS CALCULATED
        sqlQuery = u"""
            select ps1_exp_id, raDeg, decDeg from ps1_pointings where subdisks_calculated = 0
        """ % locals()

        rows = readquery(
            log=self.log,
            sqlQuery=sqlQuery,
            dbConn=self.ligo_virgo_wavesDbConn,
            quiet=False
        )
        ps1PointNum = len(rows)

        # CALCULATE ALL OF THE SUBDISKS
        inserts = []
        expIds = []
        for row in rows:
            subDiskCoordinates = self._get_subdisk_parameters(
                row["raDeg"], row["decDeg"], 1.5)
            ps1_exp_id = row["ps1_exp_id"]
            expIds.append(ps1_exp_id)
            for i, c in enumerate(subDiskCoordinates):
                insert = {
                    "raDeg": c[0],
                    "decDeg": c[1],
                    "ps1_exp_id": ps1_exp_id,
                    "circleId": i + 1
                }
                inserts.append(insert)

        # ADD SUBDISKS TO DATABASE
        if len(inserts):

            insert_list_of_dictionaries_into_database_tables(
                dbConn=self.ligo_virgo_wavesDbConn,
                log=self.log,
                dictList=inserts,
                dbTableName="ps1_pointings_subdisks",
                uniqueKeyList=["ps1_exp_id", "circleId"],
                dateModified=False,
                batchSize=2500
            )

            # UPDATE POINTINGS TABLE TO INDICATE SUBDISKS HAVE BEEN CALCULATED
            theseIds = ",".join(expIds)
            sqlQuery = u"""
                update ps1_pointings set subdisks_calculated = 1 where ps1_exp_id in (%(theseIds)s)
            """ % locals()
            writequery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
            )

        if ps1PointNum == 0:
            print "All PS1 pointings have been split into their 49 sub-disks" % locals()
        else:
            print "%(ps1PointNum)s new PS1 pointings have been split into 49 sub-disks - parameters added to the `ps1_pointings_subdisks` database table" % locals()

        # APPEND HTMIDs TO THE ps1_pointings_subdisks TABLE
        add_htm_ids_to_mysql_database_table(
            raColName="raDeg",
            declColName="decDeg",
            tableName="ps1_pointings_subdisks",
            dbConn=self.ligo_virgo_wavesDbConn,
            log=self.log,
            primaryIdColumnName="primaryId"
        )

        self.log.info(
            'completed the ``populate_ps1_subdisk_table`` method')
        return None

    def _get_subdisk_parameters(
            self,
            raDeg,
            decDeg,
            radiusDeg):
        """
        *get subdisk parameters*

        **Key Arguments:**
            - ``raDeg`` -- the central ra of the ps1 pointing
            - ``decDeg`` -- the central dec of the ps1 pointing
            - ``radiusDeg`` -- the radius of the ps1 FOV

        **Return:**
            - ``subDiskCoordinates`` -- the coordinates for 49 subdisks covering the PS1 pointing
        """
        self.log.info('starting the ``_get_subdisk_parameters`` method')

        import math
        footprintCoords = (raDeg, decDeg)

        shifts = [
            (0, 0),
            (0, math.sqrt(3.) / 2.),
            (3. / 4., math.sqrt(3.) / 4.),
            (3. / 4., -math.sqrt(3.) / 4.),
            (0, -math.sqrt(3.) / 2.),
            (-3. / 4., -math.sqrt(3.) / 4.),
            (-3. / 4., math.sqrt(3.) / 4.)
        ]

        subDiskCoordinates = []
        count = 0
        radius2 = radiusDeg / 2
        for s in shifts:
            x1 = footprintCoords[0] + s[0] * radiusDeg
            y1 = footprintCoords[1] + s[1] * radiusDeg
            for ss in shifts:
                count += 1
                x2 = x1 + ss[0] * radius2
                y2 = y1 + ss[1] * radius2
                self.log.debug("""%(count)s: %(x2)s, %(y2)s""" % locals())
                subDiskCoordinates.append((x2, y2))

        self.log.info('completed the ``_get_subdisk_parameters`` method')
        return subDiskCoordinates

    def update_ned_database_table(
            self):
        """
        *Use Sherlock & Neddy to query NED and update the catalogues database for previously unseen/stale PS1 footprint areas*

        **Return:**
            - None

        **Usage:**

            .. code-block:: python 

                # UPDATE THE NED STREAM FOR NEW PS1 FOOTPRINTS 
                from breaker import update_ps1_footprint_table
                dbUpdater = update_ps1_footprint_table(
                    log=log,
                    settings=settings
                )
                dbUpdater.update_ned_database_table()
        """
        self.log.info('starting the ``update_ned_database_table`` method')

        from sherlock.update_ned_stream import update_ned_stream

        numDisksToConesearch = 100
        rowCount = 100

        while rowCount > 0:

            sqlQuery = u"""
                select primaryId, raDeg as "ra", decDeg as "dec", htm16ID from ps1_pointings_subdisks where nedQueried = 0 limit %(numDisksToConesearch)s
            """ % locals()
            rows = readquery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
                quiet=False
            )
            rowCount = len(rows)
            ids = []
            ids[:] = [str(row["primaryId"]) for row in rows]
            ids = ",".join(ids)

            if rowCount > 0:
                print "Selecting the next %(rowCount)s subdisks areas to conesearch against NED from the `ps1_pointings_subdisks` table" % locals()
            else:
                print "NED stream is up-to-date, no queries required" % locals()

            update_ned_stream(
                log=self.log,
                cataloguesDbConn=self.cataloguesDbConn,
                settings=self.settings,
                transientsMetadataList=rows
            ).get()

            if len(ids):
                sqlQuery = u"""
                    update ps1_pointings_subdisks set nedQueried = 1 where primaryId in (%(ids)s)
                """ % locals()
                writequery(
                    log=self.log,
                    sqlQuery=sqlQuery,
                    dbConn=self.ligo_virgo_wavesDbConn,
                )

            sqlQuery = u"""
                select count(*) as count from ps1_pointings_subdisks where nedQueried = 0
            """ % locals()
            count = readquery(
                log=self.log,
                sqlQuery=sqlQuery,
                dbConn=self.ligo_virgo_wavesDbConn,
                quiet=False
            )
            count = count[0]["count"]

            if rowCount > 0:
                print "NED stream updated for %(rowCount)s PS1 pointing sub-disks (%(count)s to go)" % locals()
                print "-----\n\n"

        self.log.info('completed the ``update_ned_database_table`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method
