# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.28.0

# EMPOP's EMP files: read a small population file, count a haplotype in it for
# the lineage match probability, and write a new file for submission.
#
# Run from the deck root:  jennifer run examples/emp.j
#
# The population is made up for the example; it is not EMPOP data. A real
# file comes from fs.readString("population.emp").
use io;
use strings;

import "../src/emp.j" as emp;
import "../src/lineage.j" as lineage;

def text as string init strings.join(
    [
        "# Example population, 6 samples",
        "#! 16024-576",
        "S1\tH1c\t1\t16519C\t263G\t315.1C\t477C",
        "S2\tH5\t1\t16304C\t263G\t315.1C\t456T",
        "S3\tH5\t2\t16304C\t263G\t315.1C\t456T",
        "S4\tH2a2a1\t1\t",
        "S5\tU5a1\t1\t16192T\t16256T\t16270T\t73G\t150T\t263G\t309.1C\t315.1C",
        "S6\tJ1c\t0\t16069T\t16126C\t73G\t185A\t228A\t263G\t295T\t315.1C\t462T\t489C"
    ],
    "\n") + "\n";

def db as emp.EmpFile init emp.readEmp($text);
io.printf(
    "%s: %d entries, N = %d (count 0 excluded)\n",
    $db.comments[0],
    len($db.entries),
    emp.empSampleSize($db));
for (def e in $db.entries) {
    io.printf(
        "  line %d: %s, %s, x%d: %s\n",
        $e.line,
        $e.haplotype.id,
        $e.haplogroup,
        $e.count,
        strings.join(lineage.mtDesignations($e.haplotype), " "));
}

# The evidence, typed over the same range as the file.
def evidence as lineage.MtHaplotype init lineage.mtHaplotype(
    "evidence",
    ["16024-576"],
    "16304C 263G 315.1C 456T");
def k as int init emp.empCount($db, $evidence);
def n as int init emp.empSampleSize($db);
io.printf(
    "\nevidence seen %d times in %d: 95%% upper bound %f|prec=4, LR %f|prec=2\n",
    $k,
    $n,
    lineage.lineageMatchProbability($k, $n),
    lineage.lineageLikelihoodRatio($k, $n));

# A new file: entries are checked as readEmp checks a line.
def h as lineage.MtHaplotype init lineage.mtHaplotype(
    "NEW1",
    ["16024-16365", "73-340"],
    "73G 16189C 263G 315.1C");
def out as emp.EmpFile init emp.EmpFile{
    comments: ["New samples, HVS-I and HVS-II"],
    entries: [emp.empEntry($h, "?", 1)]
};
io.printf("\n%s", strings.replace(emp.writeEmp($out), "\t", " | "));

# A malformed line is refused with its line number.
try {
    emp.readEmp("#! 16024-576\nS1\tH\t1\t1000C\n");
} catch (e) {
    io.printf("\nrefused (%s): %s\n", $e.kind, $e.message);
}
