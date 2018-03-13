#!/usr/bin/env perl

use strict;
use warnings;

use Module::Build;

my $build = Module::Build->new(
    module_name => 'HMF::Pipeline',
    license => 'mit',
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
        'Sort::Key::Natural' => 0,
        'Template' => 0,
    },
    build_requires => {
        'Devel::Cover' => 0,
        'File::Touch' => 0,
        'Test::Cmd' => 0,
        'Test::Compile' => 0,
        'Test::Dir' => 0,
        'Test::Fatal' => 0,
        'Test::File' => 0,
        'Test::Files' => 0,
        'Test::MockModule' => 0,
        'Test::Output' => 0,
        'Task::Perl::Critic' => 0,
        'Test::Perl::Critic' => 0,
        'Test::PerlTidy' => 0,
        'Test::Prereq::Build' => '< 2.0',
        'Test::Warn' => 0,
        'Module::FromPerlVer' => 0,
    }
);

$build->create_build_script;
