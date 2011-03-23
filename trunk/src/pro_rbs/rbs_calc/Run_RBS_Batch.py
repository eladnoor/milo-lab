#!/usr/bin/python
#This file is part of the Ribosome Binding Site Calculator.

#The Ribosome Binding Site Calculator is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#The Ribosome Binding Site Calculator is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with Ribosome Binding Site Calculator.  If not, see <http://www.gnu.org/licenses/>.
#
#Copyright 2008-2009 Howard Salis
#
# - Adapted by Elad Noor, November 2009
# - Now, runs as bacth, by reading a CSV file with the sequences and start positions
# - And returns a CSV file as the result with the rest of the parameters
################################################################################

from RBS_Calculator import RBS_Calculator
import sys, math, csv

if __name__ == "__main__":

    if (len(sys.argv) < 2):
        print "Usage: " + sys.argv[0] + " [input CSV]\n"
        sys.exit(-1)
        
    input_csv = csv.reader(open(sys.argv[1], "r"))
    output_csv = csv.writer(open("result.csv", "w"))

    output_csv.writerow(('name', 'sequence', 'start position', 'expression level', 'kinetic score'))
    for row in input_csv:
        (name, seq, start) = row
        if (start == None):
            start_range = [0, len(seq)]
        else: # the new convension is to count the start position from the end of the sequence
            start_range = [int(start)-1, int(start)]
        
        #Create instance of RBS Calculator
        calcObj = RBS_Calculator(seq, start_range, name)
        calcObj.calc_dG()

        dG_total_list = calcObj.dG_total_list[:]
        start_pos_list = calcObj.start_pos_list[:]
        kinetic_score_list = calcObj.kinetic_score_list[:]

        if (len(dG_total_list) > 0):
            expr_list = []
            for dG in dG_total_list:
                expr_list.append(calcObj.K * math.exp(-dG/calcObj.RT_eff))

            for (start_pos, expr, ks) in zip(start_pos_list, expr_list, kinetic_score_list):
                output_csv.writerow((name, seq, start_pos, expr, ks))
        else:
            output_csv.writerow((name, seq, "FAIL", "FAIL", "FAIL"))            

        calcObj.cleanup()

