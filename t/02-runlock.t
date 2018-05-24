#!/usr/bin/env perl

use discipline;

use File::Temp;
use Test::Fatal;
use Test::File;
use Test::More;

use HMF::Pipeline qw(lockRun);

my $temp_dir = File::Temp->newdir();

my $lock_file = lockRun($temp_dir);
file_exists_ok($lock_file);

my $exception = exception { lockRun($temp_dir) };
like(
    $exception, qr/\QCouldn't obtain lock file (error: $!), are you *sure* there are no more jobs running? You MUST check qstat!
Jobs could have been scheduled before the error, or by external tools (GATK etc.)\E/, 'fails on already locked run'
);

done_testing();
