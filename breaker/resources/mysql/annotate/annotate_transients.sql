INSERT IGNORE INTO `tcs_gravity_event_annotations` (transient_object_id, gravity_event_id, gracedb_id) 
    select id as "transient_object_id", gravity_event_id, gracedb_id from (
        select t.id, s.earliest_mjd from tcs_latest_object_stats s, tcs_transient_objects t where s.id=t.id and t.detection_list_id != 0
    ) t, tcs_gravity_events e where t.earliest_mjd between e.mjd-10.0 and e.mjd+21.0;


DROP PROCEDURE IF EXISTS `update_atlas`;
DELIMITER //
CREATE PROCEDURE `update_atlas`()
COMMENT 'A procedure to update atlas'
BEGIN
    SELECT @exists := count(*) as count FROM information_schema.tables where table_name = 'atlas_diff_objects';
    IF @exists = 1 THEN
      INSERT IGNORE INTO `tcs_gravity_event_annotations` (transient_object_id, gravity_event_id, gracedb_id) 
          SELECT 
              id AS 'transient_object_id', gravity_event_id, gracedb_id
          FROM
              (SELECT 
                  t.id, s.earliest_mjd, t.followup_flag_date
              FROM
                  tcs_latest_object_stats s, atlas_diff_objects t
              WHERE
                  s.id = t.id
                      AND t.detection_list_id IN (1 ,2, 3, 4, 5)) t,
              tcs_gravity_events e
          WHERE
              ((TO_SECONDS(t.followup_flag_date)/(3600*24)-678941) BETWEEN e.mjd - 10.0 AND e.mjd + 21.0)
              OR (t.earliest_mjd between e.mjd-10.0 and e.mjd+21.0);
    END IF;
END//

-- CALL THE PROCEDURE
CALL `update_atlas`();
DROP PROCEDURE IF EXISTS `update_atlas`;
