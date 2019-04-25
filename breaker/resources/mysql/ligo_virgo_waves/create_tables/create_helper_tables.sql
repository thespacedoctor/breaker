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
  `prob_coverage` DOUBLE DEFAULT NULL,
  `map_name` varchar(30) DEFAULT NULL,
  `dateCreated` datetime DEFAULT NULL,
  `dateLastModified` datetime DEFAULT NULL,
  `updated` tinyint(4) DEFAULT '0',
  PRIMARY KEY (`primaryId`),
  UNIQUE KEY `skycell_id_gracedb_id` (`skycell_id`,`gracedb_id`)
) ENGINE=MyISAM AUTO_INCREMENT=0 DEFAULT CHARSET=utf8 COLLATE=utf8_swedish_ci;

CREATE TABLE IF NOT EXISTS  `atlas_exposure_gravity_event_annotations` (
  `primaryId` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `atlas_object_id` varchar(30) COLLATE utf8_swedish_ci NOT NULL,
  `gravity_event_id` varchar(10) COLLATE utf8_swedish_ci NOT NULL,
  `raDeg` DOUBLE NULL DEFAULT NULL,
  `decDeg` DOUBLE NULL DEFAULT NULL,
  `gracedb_id` varchar(10) COLLATE utf8_swedish_ci NOT NULL,
  `prob_coverage` double DEFAULT NULL,
  `map_name` varchar(30) COLLATE utf8_swedish_ci DEFAULT NULL,
  `dateCreated` datetime DEFAULT NULL,
  `dateLastModified` datetime DEFAULT NULL,
  `updated` tinyint(4) DEFAULT '0',
  PRIMARY KEY (`primaryId`),
  UNIQUE KEY `atlas_object_id_gracedb_id` (`atlas_object_id`,`gracedb_id`),
  KEY `idx_atlas_object_id` (`atlas_object_id`),
  KEY `idx_gracedb_id` (`gracedb_id`)
) ENGINE=MyISAM AUTO_INCREMENT=0 DEFAULT CHARSET=utf8 COLLATE=utf8_swedish_ci;
