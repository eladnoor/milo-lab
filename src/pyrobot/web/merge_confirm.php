<?php

$con=mysqli_connect("hldbv02","ronm","a1a1a1","tecan");
// Check connection
if (mysqli_connect_errno())
{
    echo "Failed to connect to MySQL: " . mysqli_connect_error();
}
if (!(array_key_exists("exp_condition", $_GET) && array_key_exists("merged_exp_id", $_GET)))
{
    echo "Must provide keys for exp_condition and merged_exp_id";
}

$query = "INSERT INTO tecan_readings (exp_id, plate, reading_label, row, col, time, measurement)
          SELECT   \"" . $_GET["merged_exp_id"] . "\", plate, reading_label, row, col, time, measurement FROM tecan_readings 
          WHERE    ". $_GET["exp_condition"] . " 
          ORDER BY time, plate, reading_label, row, col";

$result = mysqli_query($con, $query);
#echo $query . "</br>";

$query = "INSERT INTO tecan_experiments VALUES (\"". $_GET["merged_exp_id"]."\", \"MERGE\", \"".$_GET["description"]."\")";

$result = mysqli_query($con, $query);
#echo $query  . "</br>";

$query = "INSERT INTO tecan_plates VALUES (\"". $_GET["merged_exp_id"]."\", 0, \"".$_GET["description"]."\", NULL, NULL)";

$result = mysqli_query($con, $query);
#echo $query  . "</br>";

echo "GREAT SUCCESS!!! The experiments have been merged to the new Exp. ID: ". $_GET["merged_exp_id"]."</br>";
echo "Select other experiments to merge: <a href='exp_list.php'>Experiment list</a></br>";
?>
