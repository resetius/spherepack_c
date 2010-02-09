while (<>) {
	$_ =~ s/pow_ri/pow_di/g;
	print $_;
}
