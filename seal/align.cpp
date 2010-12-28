
/*
 * align.cpp
 * Copyright (C) Carpov Pavel   2010 <carpovpv@qsar.chem.msu.ru>
                 Baskin Igor I. 2010 <igbaskin@gmail.com>
 *
 * MCMF is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MCMF is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "seal.h"
#include "ks_seal.h"
#include <string>
#include <stdexcept>
#include <string.h>
#include <iostream>
#include <stdlib.h>

int main(int argc , char **argv)
{
    double wE = 1.0;
    double wS = 1.0;
    double alpha = 0.5;
    char * fp = NULL;
    int nprobes = 10;

    while(argc-- >0)
    {
        if(!strncmp(argv[argc],"--sdf=",6 ))
            fp = argv[argc] + 6;
        else if(!strncmp(argv[argc],"--we=", 5))
            wE = atof(argv[argc] + 5);
        else if(!strncmp(argv[argc],"--ws=", 5))
            wS = atof(argv[argc] + 5);
        else if(!strncmp(argv[argc],"--alpha=", 8))
            alpha = atof(argv[argc] + 8);
        else if(!strncmp(argv[argc], "--np=",5))
            nprobes = atoi(argv[argc] + 5);

    }

    std::cerr << "Filename: " << fp << std::endl;
    std::cerr << "Alpha: " << alpha << std::endl;
    std::cerr << "wE: " << wE << std::endl;
    std::cerr << "wS: " << wS << std::endl;
    std::cerr << "Nprobes: " << nprobes << std::endl;

    try
    {
        SEAL * seal = new KS_Seal(fp,
                                  alpha,
                                  wE,
                                  wS,
                                  nprobes);
        seal->go();

        delete seal;
    }
    catch(std::exception &exc)
    {
        fprintf(stderr, "%s\n", exc.what());
    }

    return 0;
}

