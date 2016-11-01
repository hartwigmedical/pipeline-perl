package discipline;

# for this module itself (lexical scope)
use 5.016_000;
use strict;
use warnings FATAL => 'all';

# use without importing anything
use feature ();
use strictures 2 ();
use indirect ();
use multidimensional ();
use bareword::filehandles ();
use autovivification ();
use mro ();

# for client code
# module->import is equivalent to client code writing "use module"
# module->unimport is equivalent to client code writing "no module"
# see Modern::Perl, strictures etc.: just being more specific here
sub import {
    my ($class) = @_;

    feature->import(':5.16');
    strict->import;

    # disagree that these should be allowed for the reasons given
    my @DISAGREE_NONFATAL = grep { exists $warnings::Offsets{$_} } (
        'newline',      # stat on nonexistent file with a newline in it
        'experimental', # no reason for these to be fatal
        'deprecated',   # unfortunately can't make these fatal
        'portable',     # everything worked fine here, just may not elsewhere
    );

    # import warnings
    # make all known (but not custom) fatal
    # make some not fatal
    # revert those we disagree with
    # disable some entirely
    warnings->import;
    warnings->import(FATAL => @strictures::WARNING_CATEGORIES);
    warnings->unimport(FATAL => @strictures::V2_NONFATAL);
    warnings->import(@strictures::V2_NONFATAL);
    warnings->import(FATAL => @DISAGREE_NONFATAL);
    warnings->unimport(@strictures::V2_DISABLE);

    # extras
    indirect->unimport(':fatal');
    multidimensional->unimport;
    bareword::filehandles->unimport;

    # defaults plus strict mode
    autovivification->unimport;
    autovivification->unimport('strict');

    # to consider but ignores existing user-provided handling
    # autodie->import;

    # intuitive, modern inheritence model
    mro::set_mro(scalar caller(), 'c3'); ## no critic (Subroutines::ProhibitCallsToUnexportedSubs)
    return;
}

1;
