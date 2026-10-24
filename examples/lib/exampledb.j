# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# SYNTHETIC illustrative allele frequencies for the examples - NOT real
# population data. Every number here is made up.
#
# The deck reads no published format: a reader deck converts one into a
# FrequencyDb - @deckyard/strider does that for STRidER XML, at
# https://github.com/deckyard/deck-strider-frequency-converter. The examples build one
# in code instead, so they depend on nothing but the deck itself. For anything
# you intend to rely on, get a real table from https://strider.online/.

import "../../src/forensicgenetics.j" as forensics;

/**
 * A synthetic eight-locus reference table.
 * @return {forensics.FrequencyDb} the example database
 */
export func build() {
    def db as forensics.FrequencyDb init forensics.db("Synthetic Example Population");
    $db = forensics.withLocus(
        $db,
        "D3S1358",
        1034,
        {
            "12": 0.2386,
            "13": 0.0509,
            "14": 0.0394,
            "15": 0.0825,
            "16": 0.1813,
            "17": 0.1327,
            "18": 0.2452,
            "19": 0.0294
        });
    $db = forensics.withLocus(
        $db,
        "vWA",
        1034,
        {
            "14": 0.1635,
            "15": 0.2070,
            "16": 0.0346,
            "17": 0.0541,
            "18": 0.1014,
            "19": 0.2139,
            "20": 0.2255
        });
    $db = forensics.withLocus(
        $db,
        "FGA",
        1034,
        {
            "18": 0.0145,
            "19": 0.3298,
            "20": 0.0387,
            "21": 0.1525,
            "22": 0.0175,
            "23": 0.1012,
            "24": 0.0821,
            "25": 0.2054,
            "26": 0.0583
        });
    $db = forensics.withLocus(
        $db,
        "D8S1179",
        1034,
        {
            "8": 0.2244,
            "9": 0.1253,
            "10": 0.1336,
            "11": 0.1246,
            "12": 0.0052,
            "13": 0.1157,
            "14": 0.0778,
            "15": 0.0411,
            "16": 0.1523
        });
    $db = forensics.withLocus(
        $db,
        "D21S11",
        1034,
        {
            "27": 0.0112,
            "28": 0.2348,
            "29": 0.0225,
            "30": 0.0070,
            "30.2": 0.0076,
            "31": 0.1513,
            "31.2": 0.2212,
            "32.2": 0.2301,
            "33.2": 0.1143
        });
    $db = forensics.withLocus(
        $db,
        "D18S51",
        1034,
        {
            "12": 0.1580,
            "13": 0.1604,
            "14": 0.0454,
            "15": 0.1358,
            "16": 0.1396,
            "17": 0.0835,
            "18": 0.2060,
            "19": 0.0649,
            "20": 0.0064
        });
    $db = forensics.withLocus(
        $db,
        "TH01",
        1034,
        {
            "5": 0.2764,
            "6": 0.1751,
            "7": 0.1137,
            "8": 0.3092,
            "9": 0.0382,
            "9.3": 0.0775,
            "10": 0.0099
        });
    $db = forensics.withLocus(
        $db,
        "D16S539",
        1034,
        {
            "8": 0.0541,
            "9": 0.1303,
            "10": 0.1140,
            "11": 0.4339,
            "12": 0.0196,
            "13": 0.2207,
            "14": 0.0274
        });
    return $db;
}
