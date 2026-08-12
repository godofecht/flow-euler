// Project Euler 954
// Port of the Python reference solver.
// Digit DP with mod-7 residue tracking using a hash table.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int W_arr[6] = {1, 3, 2, 6, 4, 5};
static int SHIFT_arr[6];
static int shifts_tab[7][6][6];
static int rot_tab[7][128];
static int mask_all[512];
static int mask_no0[512];
static int update_tab[512][8];
static int add_contrib_tab[6][8];
static int perm_tab[7][7];
static int res_of_idx[8] = {0, 0, 1, 2, 3, 4, 5, 6};
static int mult_of_idx[8] = {1, 1, 2, 2, 1, 1, 1, 1};

static int mod_inv7(int a) {
    for (int i = 1; i < 7; i++)
        if ((a * i) % 7 == 1) return i;
    return 0;
}

static void init_tables(void) {
    for (int c = 0; c < 6; c++) SHIFT_arr[c] = 9 * c;

    int invdiff[6][6];
    memset(invdiff, 0, sizeof(invdiff));
    for (int a = 0; a < 6; a++)
        for (int b = 0; b < 6; b++)
            if (a != b)
                invdiff[a][b] = mod_inv7((W_arr[b] - W_arr[a] + 7) % 7);

    for (int r = 1; r < 7; r++)
        for (int a = 0; a < 6; a++)
            for (int b = 0; b < 6; b++)
                if (a != b)
                    shifts_tab[r][a][b] = (r * invdiff[a][b]) % 7;

    for (int sh = 0; sh < 7; sh++)
        for (int m = 0; m < 128; m++)
            rot_tab[sh][m] = (sh == 0) ? m : ((m << sh) | (m >> (7 - sh))) & 0x7F;

    for (int bits = 0; bits < 512; bits++) {
        int mask = bits & 0x7F;
        int has7 = (bits >> 8) & 1;
        mask_all[bits] = mask;
        int m2 = mask & ~1;
        if (has7) m2 |= 1;
        mask_no0[bits] = m2;
    }

    int choices_by_idx[8][3] = {
        {0,1,0},{0,0,1},{1,0,0},{2,0,0},
        {3,0,0},{4,0,0},{5,0,0},{6,0,0}
    };
    for (int oldbits = 0; oldbits < 512; oldbits++) {
        int mask = oldbits & 0x7F;
        int has0 = (oldbits >> 7) & 1;
        int has7 = (oldbits >> 8) & 1;
        for (int idx = 0; idx < 8; idx++) {
            int res = choices_by_idx[idx][0];
            int add0 = choices_by_idx[idx][1];
            int add7 = choices_by_idx[idx][2];
            int newmask = mask | (1 << res);
            int newhas0 = has0 | add0;
            int newhas7 = has7 | add7;
            update_tab[oldbits][idx] = newmask | (newhas0 << 7) | (newhas7 << 8);
        }
    }

    for (int c = 0; c < 6; c++)
        for (int idx = 0; idx < 8; idx++)
            add_contrib_tab[c][idx] = (res_of_idx[idx] * W_arr[c]) % 7;

    for (int add = 0; add < 7; add++)
        for (int i = 0; i < 7; i++)
            perm_tab[add][i] = (i + add) % 7;
}

// Hash table with open addressing, linear probing, and slot list
typedef struct {
    int capacity;
    int mask;
    int count;
    char *occupied;
    long long *states;
    int *counts;  // capacity * 7
    int *slot_list;
    int slot_list_size;
} HT;

static unsigned int ht_hash(long long state, int mask) {
    unsigned long long x = (unsigned long long)state;
    x ^= x >> 33;
    x *= 0xFF51AFD7ED558CCDULL;
    x ^= x >> 33;
    x *= 0xC4CEB9FE1A85EC53ULL;
    x ^= x >> 33;
    return (unsigned int)(x & (unsigned int)mask);
}

static void ht_init(HT *ht, int capacity) {
    ht->capacity = capacity;
    ht->mask = capacity - 1;
    ht->count = 0;
    ht->occupied = (char *)calloc(capacity, 1);
    ht->states = (long long *)malloc(capacity * sizeof(long long));
    ht->counts = (int *)malloc((size_t)capacity * 7 * sizeof(int));
    ht->slot_list = (int *)malloc(capacity * sizeof(int));
    ht->slot_list_size = 0;
}

static void ht_clear(HT *ht) {
    for (int i = 0; i < ht->slot_list_size; i++) {
        ht->occupied[ht->slot_list[i]] = 0;
    }
    ht->count = 0;
    ht->slot_list_size = 0;
}

static void ht_free(HT *ht) {
    free(ht->occupied);
    free(ht->states);
    free(ht->counts);
    free(ht->slot_list);
}

static void ht_resize(HT *ht, int new_cap) {
    HT old = *ht;
    ht->capacity = new_cap;
    ht->mask = new_cap - 1;
    ht->count = 0;
    ht->slot_list_size = 0;
    ht->occupied = (char *)calloc(new_cap, 1);
    ht->states = (long long *)malloc(new_cap * sizeof(long long));
    ht->counts = (int *)malloc((size_t)new_cap * 7 * sizeof(int));
    ht->slot_list = (int *)malloc(new_cap * sizeof(int));

    for (int si = 0; si < old.slot_list_size; si++) {
        int slot = old.slot_list[si];
        long long state = old.states[slot];
        unsigned int h = ht_hash(state, ht->mask);
        while (ht->occupied[h]) {
            h = (h + 1) & ht->mask;
        }
        ht->occupied[h] = 1;
        ht->states[h] = state;
        int *c = &ht->counts[(size_t)h * 7];
        int *oc = &old.counts[(size_t)slot * 7];
        for (int j = 0; j < 7; j++) c[j] = oc[j];
        ht->slot_list[ht->slot_list_size++] = h;
        ht->count++;
    }

    free(old.occupied);
    free(old.states);
    free(old.counts);
    free(old.slot_list);
}

static inline void ht_insert(HT *ht, long long state, const int *cnts) {
    if (ht->count * 2 >= ht->capacity) {
        ht_resize(ht, ht->capacity * 2);
    }
    unsigned int h = ht_hash(state, ht->mask);
    while (ht->occupied[h]) {
        if (ht->states[h] == state) {
            int *c = &ht->counts[(size_t)h * 7];
            c[0] += cnts[0]; c[1] += cnts[1]; c[2] += cnts[2];
            c[3] += cnts[3]; c[4] += cnts[4]; c[5] += cnts[5];
            c[6] += cnts[6];
            return;
        }
        h = (h + 1) & ht->mask;
    }
    ht->occupied[h] = 1;
    ht->states[h] = state;
    int *c = &ht->counts[(size_t)h * 7];
    c[0] = cnts[0]; c[1] = cnts[1]; c[2] = cnts[2];
    c[3] = cnts[3]; c[4] = cnts[4]; c[5] = cnts[5];
    c[6] = cnts[6];
    ht->slot_list[ht->slot_list_size++] = (int)h;
    ht->count++;
}

static void advance_dp(HT *old_ht, HT *new_ht, int pos, int target_r, int is_MSD) {
    int c = pos % 6;
    int shiftc = SHIFT_arr[c];
    int *mask_func = is_MSD ? mask_no0 : mask_all;
    int choices_start = is_MSD ? 1 : 0;
    int sh_row[6];
    for (int a = 0; a < 6; a++) sh_row[a] = shifts_tab[target_r][a][c];

    ht_clear(new_ht);

    for (int si = 0; si < old_ht->slot_list_size; si++) {
        int i = old_ht->slot_list[si];
        long long state = old_ht->states[i];
        const int *cnts = &old_ht->counts[(size_t)i * 7];

        int forb = 0;
        for (int a = 0; a < 6; a++) {
            if (a == c) continue;
            int bitsa = (int)((state >> SHIFT_arr[a]) & 0x1FF);
            int mask_use = mask_func[bitsa];
            if (mask_use) {
                forb |= rot_tab[sh_row[a]][mask_use];
            }
        }

        int bitsc = (int)((state >> shiftc) & 0x1FF);

        for (int idx = choices_start; idx < 8; idx++) {
            int res = res_of_idx[idx];
            if (forb & (1 << res)) continue;
            int newbitsc = update_tab[bitsc][idx];
            long long newstate = state ^ ((long long)(bitsc ^ newbitsc) << shiftc);
            int add = add_contrib_tab[c][idx];
            int mult = mult_of_idx[idx];

            int newcnts[7];
            const int *p = perm_tab[add];
            if (mult == 1) {
                newcnts[p[0]] = cnts[0];
                newcnts[p[1]] = cnts[1];
                newcnts[p[2]] = cnts[2];
                newcnts[p[3]] = cnts[3];
                newcnts[p[4]] = cnts[4];
                newcnts[p[5]] = cnts[5];
                newcnts[p[6]] = cnts[6];
            } else {
                newcnts[p[0]] = cnts[0] * 2;
                newcnts[p[1]] = cnts[1] * 2;
                newcnts[p[2]] = cnts[2] * 2;
                newcnts[p[3]] = cnts[3] * 2;
                newcnts[p[4]] = cnts[4] * 2;
                newcnts[p[5]] = cnts[5] * 2;
                newcnts[p[6]] = cnts[6] * 2;
            }

            ht_insert(new_ht, newstate, newcnts);
        }
    }
}

static long long count_len_res(int L, int target_r, HT *ht_a, HT *ht_b) {
    ht_clear(ht_a);
    int init_cnts[7] = {1, 0, 0, 0, 0, 0, 0};
    ht_insert(ht_a, 0, init_cnts);

    for (int pos = 0; pos < L; pos++) {
        advance_dp(ht_a, ht_b, pos, target_r, pos == L - 1);
        HT tmp = *ht_a; *ht_a = *ht_b; *ht_b = tmp;
    }

    long long total = 0;
    for (int si = 0; si < ht_a->slot_list_size; si++) {
        int i = ht_a->slot_list[si];
        total += ht_a->counts[(size_t)i * 7 + target_r];
    }
    return total;
}

long long p954_native(void) {
    init_tables();

    HT ht_a, ht_b;
    ht_init(&ht_a, 1024);
    ht_init(&ht_b, 1024);

    long long total = 0;
    for (int L = 1; L <= 13; L++) {
        for (int r = 1; r <= 6; r++) {
            total += count_len_res(L, r, &ht_a, &ht_b);
        }
    }

    ht_free(&ht_a);
    ht_free(&ht_b);

    return total;
}
