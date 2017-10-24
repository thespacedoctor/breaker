CREATE TABLE IF NOT EXISTS `gravity_events` (
  `gravity_event_id` varchar(10) COLLATE utf8_unicode_ci NOT NULL,
  `gracedb_id` varchar(10) COLLATE utf8_unicode_ci NOT NULL,
  `mjd` DOUBLE DEFAULT NULL,
  `dateLastModified` datetime DEFAULT NULL,
  `dateCreated` datetime DEFAULT NULL,
  `updated` tinyint(4) DEFAULT '0',
  PRIMARY KEY (`gravity_event_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_swedish_ci;


CREATE TABLE IF NOT EXISTS `ps1_skycell_gravity_event_annotations` (
  `primaryId` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `skycell_id` varchar(30) NOT NULL,
  `gravity_event_id` varchar(30) NOT NULL,
  `gracedb_id` varchar(10) NOT NULL,
  `enclosing_contour` int(11) DEFAULT NULL,
  `map_name` varchar(30) DEFAULT NULL,
  `dateCreated` datetime DEFAULT NULL,
  `dateLastModified` datetime DEFAULT NULL,
  `updated` tinyint(4) DEFAULT '0',
  PRIMARY KEY (`primaryId`),
  UNIQUE KEY `skycell_id_gracedb_id` (`skycell_id`,`gracedb_id`)
) ENGINE=MyISAM AUTO_INCREMENT=0 DEFAULT CHARSET=utf8 COLLATE=utf8_swedish_ci;


INSERT IGNORE INTO `ps1_skycell_gravity_event_annotations` (skycell_id, gravity_event_id, gracedb_id) 
SELECT 
    distinct skycell_id, gravity_event_id, gracedb_id
FROM
    (SELECT 
        skycell_id, mjd
    FROM
        ps1_warp_stack_diff_skycells UNION ALL SELECT 
        skycell_id, mjd
    FROM
        ps1_stack_stack_diff_skycells) t,
    gravity_events e
WHERE
    t.mjd BETWEEN e.mjd - 21.0 AND e.mjd + 35.0;
