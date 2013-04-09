<?php
#date_default_timezone_set('Etc/GMT');

$con=mysqli_connect("hldbv02","ronm","a1a1a1","tecan");
// Check connection
if (mysqli_connect_errno())
{
    echo "Failed to connect to MySQL: " . mysqli_connect_error();
}

if (!array_key_exists("exp_id", $_GET))
{
    $result = mysqli_query($con,"SELECT exp_id, serial_number, description FROM tecan_experiments ORDER BY exp_id DESC");

    echo "<form name=\"input\" action=\"merge.php\" method=\"get\">";
    echo "Merged exp_id: <input type=\"text\" name=\"merged_exp_id\">";
    echo "<input type=\"submit\" value=\"Submit\"></br>";

    echo "<table border='1'>\n";
    echo "<tr>\n";
    echo "<th></th>\n";
    echo "<th>Exp. ID</th>\n";
    echo "<th>Description</th>\n";
    echo "</tr>\n";

    while($row = mysqli_fetch_array($result))
    {
        $href = "\"?exp_id=" . $row['exp_id'] . "&row=0&col=0&reading_label=OD600&plate=0\"";
        echo "<tr>";
        echo "<td><input type=\"checkbox\" name=\"exp_id[]\" value=\"". $row['exp_id'] ."\"></td>";
        echo "<td><a href=" . $href . ">" . $row['exp_id'] . "</a></td>";
        echo "<td>" . $row['description'] . "</td>\n";
        echo "</tr>";
    }
    echo "</table>";
}
else
{
    if (!(array_key_exists("plate", $_GET) && array_key_exists("row", $_GET) && array_key_exists("col", $_GET) && array_key_exists("reading_label", $_GET)))
    {
        echo "ERROR: must provide plate, reading_label, row and col</br>";
    }
    else
    {
        $start_time = strtotime($_GET["exp_id"]);
        $query = "SELECT time, measurement FROM tecan_readings 
                  WHERE exp_id=\"" . $_GET["exp_id"] . "\" 
                  AND reading_label=\"" . $_GET["reading_label"] . "\" 
                  AND plate=" . $_GET["plate"] . " 
                  AND row=" . $_GET["row"] . " 
                  AND col=" . $_GET["col"];

        $result = mysqli_query($con,$query);
        echo "<table border='1'>";
        echo "<tr>";
        echo "<th>Time (minutes)</th>";
        echo "<th>Measurement</th>";
        echo "</tr>\n";

        while($row = mysqli_fetch_array($result))
        {
            echo "<tr>";
            echo "<td>" . round(($row['time'] - $start_time)/60) . "</td>";
            echo "<td>" . $row['measurement'] . "</td>";
            echo "</tr>\n";
        }
        echo "</table>";
    }
}
mysqli_close($con);
?>
