while (<>) {
	$_ =~ s/([ (	{,.*]+)real([ )},.*	]+)/\1doublereal\2/g;
	print $_;
}
