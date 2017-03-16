CREATE TABLE IF NOT EXISTS `tcs_gravity_events` (
  `gravity_event_id` varchar(10) COLLATE utf8_unicode_ci NOT NULL,
  `gracedb_id` varchar(10) COLLATE utf8_unicode_ci NOT NULL,
  `mjd` DOUBLE DEFAULT NULL,
  `dateLastModified` datetime DEFAULT NULL,
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
  `dateLastModified` datetime DEFAULT NULL,
  `updated` tinyint(4) DEFAULT '0',
  PRIMARY KEY (`primaryId`),
  UNIQUE KEY `transient_gracedb` (`transient_object_id`,`gracedb_id`)
) ENGINE=MyISAM AUTO_INCREMENT=0 DEFAULT CHARSET=utf8 COLLATE=utf8_swedish_ci;


INSERT IGNORE INTO `tcs_gravity_event_annotations` (transient_object_id, gravity_event_id, gracedb_id) 
    select id as "transient_object_id", gravity_event_id, gracedb_id from (
        select t.id, s.earliest_mjd from tcs_latest_object_stats s, tcs_transient_objects t where s.id=t.id and t.detection_list_id != 0
    ) t, tcs_gravity_events e where t.earliest_mjd between e.mjd-21.0 and e.mjd+35.0;


DROP PROCEDURE IF EXISTS `update_atlas`;
DELIMITER //
CREATE PROCEDURE `update_atlas`()
COMMENT 'A procedure to update atlas'
BEGIN
    SELECT @exists := count(*) as count FROM information_schema.tables where table_name = 'atlas_diff_objects';
    IF @exists = 1 THEN
      INSERT IGNORE INTO `tcs_gravity_event_annotations` (transient_object_id, gravity_event_id, gracedb_id) 
          select id as "transient_object_id", gravity_event_id, gracedb_id from (
              select t.id, s.earliest_mjd from tcs_latest_object_stats s, atlas_diff_objects t where s.id=t.id and t.detection_list_id in (1,2,3,4,5)
          ) t, tcs_gravity_events e where t.earliest_mjd between e.mjd-21.0 and e.mjd+35.0;
    END IF;
END//

-- CALL THE PROCEDURE
CALL `update_atlas`();
DROP PROCEDURE IF EXISTS `update_atlas`;
