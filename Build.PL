#!/usr/bin/env perl

use strict;
use warnings;

use Module::Build;

my $build = Module::Build->new(
    module_name => 'HMF::Pipeline',
    license  => 'mit',
    requires => {
        'perl' => '5.16.0',
        'strictures' => '2.0',
        'indirect' => 0,
        'multidimensional' => 0,
        'bareword::filehandles' => 0,
        'autovivification' => 0,
        'File::Copy::Recursive' => 0,
        'File::Find::Rule' => 0,
        'FindBin::libs' => 0,
        'JSON' => 0,
        'List::MoreUtils' => 0,
        'Template' => 0,
    },
    build_requires => {
        'Test::Compile' => 0,
        'Task::Perl::Critic' => 0,
        'Test::Perl::Critic' => 0,
        'Test::PerlTidy' => 0,
        'Test::Prereq::Build' => '< 2.0',
    }
);

$build->create_build_script;
