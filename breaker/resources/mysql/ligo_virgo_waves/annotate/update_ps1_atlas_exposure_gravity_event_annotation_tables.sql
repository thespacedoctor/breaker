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
    t.mjd BETWEEN e.mjd - 10.0 AND e.mjd + 21.0;


INSERT IGNORE INTO `atlas_exposure_gravity_event_annotations` (atlas_object_id, raDEg, decDeg, gravity_event_id, gracedb_id) 
SELECT 
    distinct atlas_object_id, raDeg, decDeg, gravity_event_id, gracedb_id
FROM
    (SELECT 
        atlas_object_id, mjd, raDeg, decDeg
    FROM
        atlas_pointings) t,
    gravity_events e
WHERE
    t.mjd BETWEEN e.mjd - 10.0 AND e.mjd + 21.0;
