#include <stdio.h>
#include <string.h>
static void append_first(const char *s, char *out, int *olen){
    int word_start=1;
    for(int i=0;s[i];i++){
        if(s[i]==' ') word_start=1;
        else if(word_start){ out[(*olen)++]=s[i]; word_start=0; }
    }
}
void pe_solve_print(void){
    char out[64]; int olen=0;
    append_first("affine plane", out, &olen);
    append_first("radically integral local field", out, &olen);
    append_first("open oriented line section", out, &olen);
    append_first("jacobian", out, &olen);
    append_first("orthogonal kernel embedding", out, &olen);
    out[olen]=0;
    printf("%s\n", out);
}
