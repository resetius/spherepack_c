while (<>) {
	$_ =~ s/( [0-9]+\.[0-9]*)f/\1/g;
	$_ =~ s/( [0-9]*\.[0-9]+)f/\1/g;
	$_ =~ s/( [0-9]+)f/\1/g;
	print $_;
}
