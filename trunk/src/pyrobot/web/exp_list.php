<?php
$con=mysqli_connect("hldbv02","ronm","a1a1a1","tecan");
// Check connection
if (mysqli_connect_errno())
  {
  echo "Failed to connect to MySQL: " . mysqli_connect_error();
  }

$result = mysqli_query($con,"SELECT exp_id, serial_number, description FROM tecan_experiments ORDER BY exp_id DESC");

echo "<table border='1'>
<tr>
<th>Firstname</th>
<th>Lastname</th>
</tr>";

while($row = mysqli_fetch_array($result))
  {
  echo "<tr>";
  echo "<td>" . $row['exp_id'] . "</td>";
  echo "<td>" . $row['description'] . "</td>";
  echo "</tr>";
  }
echo "</table>";

mysqli_close($con);
?>
