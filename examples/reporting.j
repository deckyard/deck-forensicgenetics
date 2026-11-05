# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

# Reporting a single-source match: the plug-in LR, the size-bias correction,
# the FSR-GUI-0012 Bayesian estimate, a 99% posterior lower bound (allowances
# for the database's sampling error), and the ENFSI (2015) verbal equivalent
# of each.
#
# Run from the deck root:  jennifer run examples/reporting.j
#
# The frequency table is synthetic; see examples/lib/exampledb.j.
use io;

import "../src/forensicgenetics.j" as forensics;
import "../src/reporting.j" as reporting;
import "./lib/exampledb.j" as exampledb;

def db as forensics.FrequencyDb init exampledb.build();
def theta as float init forensics.THETA_GENERAL;

def suspect as forensics.Profile init forensics.profile(
    "suspect",
    [
        forensics.genotype("D3S1358", "15", "16"),
        forensics.genotype("vWA", "17", "17"),
        forensics.genotype("FGA", "21", "24"),
        forensics.genotype("D8S1179", "12", "14"),
        forensics.genotype("TH01", "6", "9.3")
    ]);

func show(name as string, lr as float) {
    def v as reporting.VerbalEquivalent init reporting.verbalEquivalent($lr);
    io.printf("%s|pad=26 LR %f|prec=0|pad=10  %s\n", $name, $lr, $v.label);
}

show("plug-in, 5/2N floor", forensics.matchLikelihoodRatio($db, $suspect, $theta));
show("size-bias corrected", reporting.sizeBiasLikelihoodRatio($db, $suspect, $theta));
show(
    "Bayesian, 1/(k+2) prior",
    reporting.bayesianLikelihoodRatio($db, $suspect, $theta, reporting.PRIOR_ONE_OVER_K_PLUS_TWO));
show(
    "99% lower bound, 1/k prior",
    reporting.likelihoodRatioLowerBound(
        $db,
        $suspect,
        $theta,
        reporting.PRIOR_ONE_OVER_K,
        0.99,
        2000,
        1));

def q as list of float init reporting.log10RmpQuantiles(
    $db,
    $suspect,
    $theta,
    reporting.PRIOR_UNIFORM,
    [0.01, 0.5, 0.99],
    2000,
    1);
io.printf(
    "\nposterior log10 RMP, uniform prior: 1%% %f|prec=3  50%% %f|prec=3  99%% %f|prec=3\n",
    $q[0],
    $q[1],
    $q[2]);

def words as reporting.VerbalEquivalent init reporting.verbalEquivalent(
    reporting.sizeBiasLikelihoodRatio($db, $suspect, $theta));
io.printf("\n%s\n", $words.statement);
