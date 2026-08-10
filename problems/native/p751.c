#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void sequence(double theta, int length, char *out) {
    double b1 = theta;
    int pos = 0;
    pos += sprintf(out + pos, "%d.", (int)floor(b1));
    while (pos < length) {
        int floorb1 = (int)floor(b1);
        double bn = floorb1 * (b1 - floorb1 + 1.0);
        char buf[32];
        int n = sprintf(buf, "%d", (int)floor(bn));
        for (int i = 0; i < n && pos < length; i++) out[pos++] = buf[i];
        b1 = bn;
    }
    out[length] = 0;
}

void pe_solve_print(void) {
    char curr[64] = "2.";
    int length = 26;
    while ((int)strlen(curr) < length) {
        int found = -1;
        for (int x = 0; x < 10; x++) {
            char temp[64];
            sprintf(temp, "%s%d", curr, x);
            char seq[64];
            sequence(strtod(temp, NULL), (int)strlen(temp), seq);
            if (strcmp(temp, seq) == 0) { found = x; break; }
        }
        if (found < 0) { printf("fail\n"); return; }
        int L = (int)strlen(curr);
        curr[L] = (char)('0' + found);
        curr[L + 1] = 0;
    }
    curr[length] = 0;
    printf("%s\n", curr);
}
