/*
 * Project Euler 766 - Sliding Block Puzzle
 *
 * BFS over reachable configurations of a 6x5 board with 14 pieces of 6 types.
 * Pieces slide any positive number of grid units orthogonally, without rotation.
 * Pieces of the same shape are indistinguishable, so a configuration is
 * identified by the sorted anchor positions within each type.
 *
 * State encoding: 14 pieces x 5 bits = 70 bits, stored in unsigned __int128.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;
typedef uint8_t u8;

#define W 6
#define H 5
#define NCELLS (W * H)   /* 30 */
#define NTYPES 6

static const int deltas[4] = { -W, W, -1, 1 };  /* up, down, left, right */

typedef struct {
    int k;                /* number of pieces of this type */
    int noffs;            /* number of (dx,dy) offsets */
    int offs[4][2];       /* offsets */
    u32  mask_at[NCELLS]; /* occupancy bitmask for each anchor pos */
    int  limits[NCELLS][4]; /* max steps up/down/left/right */
    int  valid[NCELLS];   /* whether anchor pos is in-bounds */
    int  shift0;          /* bit offset of this type's segment */
    int  seg_bits;        /* width of this type's segment in bits */
    u128 seg_mask;        /* bitmask isolating this type's segment */
} PieceType;

static PieceType g_types[NTYPES];

static void init_types(void) {
    /* offsets, count, noffs for each type */
    static const int offs_tbl[NTYPES][4][2] = {
        {{0,0},{0,1},{1,0},{0,0}},  /* Red L triomino (3) */
        {{0,1},{1,0},{1,1},{0,0}},  /* Green L triomino (3) */
        {{0,0},{0,1},{0,0},{0,0}},  /* Yellow vertical domino (2) */
        {{0,0},{0,0},{0,0},{0,0}},  /* Magenta 1x1 (1) */
        {{0,0},{1,0},{0,1},{1,1}},  /* Blue 2x2 (4) */
        {{0,0},{1,0},{0,0},{0,0}},  /* Cyan horizontal domino (2) */
    };
    static const int noffs_tbl[NTYPES] = { 3, 3, 2, 1, 4, 2 };
    static const int count_tbl[NTYPES] = { 2, 2, 2, 6, 1, 1 };

    int start = 0;
    for (int t = 0; t < NTYPES; t++) {
        PieceType *pt = &g_types[t];
        pt->k = count_tbl[t];
        pt->noffs = noffs_tbl[t];
        for (int i = 0; i < pt->noffs; i++) {
            pt->offs[i][0] = offs_tbl[t][i][0];
            pt->offs[i][1] = offs_tbl[t][i][1];
        }
        for (int pos = 0; pos < NCELLS; pos++) {
            int x = pos % W, y = pos / W;
            u32 m = 0;
            int ok = 1;
            for (int i = 0; i < pt->noffs; i++) {
                int xx = x + pt->offs[i][0];
                int yy = y + pt->offs[i][1];
                if (xx < 0 || xx >= W || yy < 0 || yy >= H) { ok = 0; break; }
                m |= 1U << (yy * W + xx);
            }
            pt->valid[pos] = ok;
            if (!ok) continue;
            pt->mask_at[pos] = m;
            int mu = W, md = W, ml = W, mr = W;
            for (int i = 0; i < pt->noffs; i++) {
                int dy = pt->offs[i][1], dx = pt->offs[i][0];
                if (y + dy < mu) mu = y + dy;
                if ((H-1)-(y+dy) < md) md = (H-1)-(y+dy);
                if (x + dx < ml) ml = x + dx;
                if ((W-1)-(x+dx) < mr) mr = (W-1)-(x+dx);
            }
            pt->limits[pos][0] = mu;
            pt->limits[pos][1] = md;
            pt->limits[pos][2] = ml;
            pt->limits[pos][3] = mr;
        }
        pt->shift0 = start * 5;
        pt->seg_bits = pt->k * 5;
        u128 mask = 0;
        for (int i = 0; i < pt->seg_bits; i++)
            mask |= ((u128)1) << (pt->shift0 + i);
        pt->seg_mask = mask;
        start += pt->k;
    }
}

/* ---- hash set for u128 keys (open addressing, linear probing) ---- */

static inline u64 mix64(u64 h) {
    h ^= h >> 33;
    h *= 0xFF51AFD7ED558CCDULL;
    h ^= h >> 33;
    h *= 0xC4CEB9FE1A85EC53ULL;
    h ^= h >> 33;
    return h;
}

static inline u64 hash128(u128 x) {
    u64 lo = (u64)x;
    u64 hi = (u64)(x >> 64);
    return mix64(lo ^ (hi + 0x9E3779B97F4A7C15ULL));
}

#define HCAP (1u << 23)   /* 8M slots, load factor ~0.33 for 2.6M entries */

long long p766_native(void) {
    init_types();

    /* initial positions per type (already sorted) */
    static const int init_pos[NTYPES][6] = {
        {1, 4},
        {2, 22},
        {11, 16},
        {12, 13, 18, 19, 24, 25},
        {14},
        {26},
    };

    u128 state0 = 0;
    for (int t = 0; t < NTYPES; t++) {
        u128 seg = 0;
        for (int i = 0; i < g_types[t].k; i++)
            seg |= ((u128)init_pos[t][i]) << (5 * i);
        state0 |= seg << g_types[t].shift0;
    }

    /* allocate hash set */
    u128 *keys = calloc(HCAP, sizeof(u128));
    u8 *occ_flag = calloc(HCAP, 1);
    u64 hmask = HCAP - 1;

    /* allocate BFS queue */
    u128 *queue = malloc(HCAP * sizeof(u128));
    u64 qhead = 0, qtail = 0;

    /* insert initial state */
    {
        u64 i = hash128(state0) & hmask;
        while (occ_flag[i]) i = (i + 1) & hmask;
        keys[i] = state0;
        occ_flag[i] = 1;
        queue[qtail++] = state0;
    }

    u64 count = 1;

    while (qhead < qtail) {
        u128 s = queue[qhead++];

        /* decode positions and compute total occupancy */
        u32 occ = 0;
        int decoded[NTYPES][6];
        for (int t = 0; t < NTYPES; t++) {
            u128 seg = (s & g_types[t].seg_mask) >> g_types[t].shift0;
            for (int i = 0; i < g_types[t].k; i++) {
                decoded[t][i] = (int)((seg >> (5 * i)) & 0x1f);
                occ |= g_types[t].mask_at[decoded[t][i]];
            }
        }

        /* try moving each piece in each direction */
        for (int t = 0; t < NTYPES; t++) {
            int k = g_types[t].k;
            for (int j = 0; j < k; j++) {
                int pos = decoded[t][j];
                u32 m_old = g_types[t].mask_at[pos];
                u32 occ_wo = occ ^ m_old;

                for (int d = 0; d < 4; d++) {
                    int limit = g_types[t].limits[pos][d];
                    if (limit <= 0) continue;
                    for (int step = 1; step <= limit; step++) {
                        int new_pos = pos + deltas[d] * step;
                        u32 new_mask = g_types[t].mask_at[new_pos];
                        if (new_mask & occ_wo) break;

                        /* rebuild sorted position list for this type */
                        int pl[6];
                        for (int i = 0; i < k; i++) pl[i] = decoded[t][i];
                        pl[j] = new_pos;
                        int idx = j;
                        while (idx > 0 && pl[idx] < pl[idx-1]) {
                            int tmp = pl[idx]; pl[idx] = pl[idx-1]; pl[idx-1] = tmp; idx--;
                        }
                        while (idx < k-1 && pl[idx] > pl[idx+1]) {
                            int tmp = pl[idx]; pl[idx] = pl[idx+1]; pl[idx+1] = tmp; idx++;
                        }

                        u128 new_seg = 0;
                        for (int i = 0; i < k; i++)
                            new_seg |= ((u128)pl[i]) << (5 * i);
                        u128 new_state = (s & ~g_types[t].seg_mask) | (new_seg << g_types[t].shift0);

                        u64 i = hash128(new_state) & hmask;
                        int found = 0;
                        while (occ_flag[i]) {
                            if (keys[i] == new_state) { found = 1; break; }
                            i = (i + 1) & hmask;
                        }
                        if (!found) {
                            keys[i] = new_state;
                            occ_flag[i] = 1;
                            count++;
                            queue[qtail++] = new_state;
                        }
                    }
                }
            }
        }
    }

    free(keys);
    free(occ_flag);
    free(queue);
    return (long long)count;
}
