#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <openbabel/fingerprint.h>

#include <iostream>
#include <vector>

using namespace OpenBabel;
using namespace std;


char* itoa(unsigned int value, char* result, int base)
{
    // check that the base if valid
    if (base < 2 || base > 36)
    {
        *result = '\0';
        return result;
    }
    char* ptr = result, *ptr1 = result, tmp_char;
    int tmp_value;

    do {
        tmp_value = value;
        value /= base;
        *ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
    } while ( value );

    // Apply negative sign
    if (tmp_value < 0) *ptr++ = '-';
    *ptr-- = '\0';
    while (ptr1 < ptr) {
        tmp_char = *ptr;
        *ptr--= *ptr1;
        *ptr1++ = tmp_char;
    }
    return result;
}


void toSvm(unsigned int r , int n)
{

    char buffer[1024];
    itoa(r, buffer, 2);

    int i;

//    printf("%s\n", buffer);
//    printf("Buffer: %d %d\n", strlen(buffer), sizeof(int));


    for ( i=0; i<32-strlen(buffer); i++)
        printf("%c ",  '0');


    for (i=0; i< strlen(buffer); ++i)
        printf("%c ", buffer[i]);



}

int main()
{

    OBConversion conv(&cin,&cout);
    conv.SetInFormat("SDF");

    OBMol mol;
    OBFingerprint *ob  = OBFingerprint::FindFingerprint("FP2");

    for (int i =0; i< 1024; i++)
        printf("_%d ", i+1);

    int Q = 0;
    printf("\n");

    while (conv.Read(&mol))
    {
        char buffer[1024];
        char zeros[1024];

        vector<unsigned int> fp;
        int N =1;

        ob->GetFingerprint(&mol,fp, 1024);


        for (int i=0; i< fp.size(); ++i, N+=32)
          toSvm(fp[i], N);
	printf("\n");

    }
    printf("$$$$ fp2 0 1023");
    return 0;

}
