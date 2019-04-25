CREATE TABLE IF NOT EXISTS `tcs_gravity_events` (
  `gravity_event_id` varchar(10) COLLATE utf8_unicode_ci NOT NULL,
  `gracedb_id` varchar(10) COLLATE utf8_unicode_ci NOT NULL,
  `mjd` DOUBLE DEFAULT NULL,
  `dateLastModified` datetime DEFAULT NULL,
  `dateCreated` datetime DEFAULT NULL,
  `updated` tinyint(4) DEFAULT '0',
  PRIMARY KEY (`gravity_event_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_swedish_ci;


CREATE TABLE IF NOT EXISTS `tcs_gravity_event_annotations` (
  `primaryId` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `transient_object_id` bigint(20) unsigned NOT NULL,
  `gravity_event_id` varchar(10) NOT NULL,
  `gracedb_id` varchar(10) NOT NULL,
  `enclosing_contour` int(11) DEFAULT NULL,
  `map_name` varchar(30) DEFAULT NULL,
  `dateCreated` datetime DEFAULT NULL,
  `dateLastModified` datetime DEFAULT NULL,
  `updated` tinyint(4) DEFAULT '0',
  PRIMARY KEY (`primaryId`),
  UNIQUE KEY `transient_gracedb` (`transient_object_id`,`gracedb_id`)
) ENGINE=MyISAM AUTO_INCREMENT=0 DEFAULT CHARSET=utf8 COLLATE=utf8_swedish_ci;


