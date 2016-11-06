#!/usr/bin/env perl

use strict;
use warnings;

use File::Temp;
use Test::Fatal;
use Test::File;
use Test::More;

use HMF::Pipeline qw(lockRun);


my $temp_dir = File::Temp->newdir();

my $lock_file = lockRun($temp_dir);
file_exists_ok($lock_file);

my $exception = exception { lockRun($temp_dir) };
like($exception, qr/\QCouldn't obtain lock file, are you *sure* there are no more jobs running?\E/, 'fails on already locked run');

done_testing();
