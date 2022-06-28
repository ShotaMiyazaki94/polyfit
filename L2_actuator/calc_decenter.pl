#! /usr/bin/perl


use Math::Trig;
() = @ARGV;

$Rx = 0;
$Ry =	0;
for($i=0;$i<@ARGV;$i++){
	$Rx = $ARGV[$i+1] if $ARGV[$i] eq "-Rx";	#	[as]
	$Ry = $ARGV[$i+1] if $ARGV[$i] eq "-Ry";	#	[as]
}
$distance_CNR_Spider = 500;	# [mm]

$shiftx = $distance_CNR_Spider*$Ry/3600.0/180.0*pi;
$shifty = $distance_CNR_Spider*$Rx/3600.0/180.0*pi;
printf "move shiftx = %.3f [mm], thus screw %.3f [spin] %.3f [deg]\n",$shiftx,$shiftx/2.5,$shiftx/2.5*360;
printf "move shifty = %.3f [mm], thus screw %.3f [spin] %.3f [deg]\n",$shifty,$shifty/2.5,$shifty/2.5*360;

exit(1);
