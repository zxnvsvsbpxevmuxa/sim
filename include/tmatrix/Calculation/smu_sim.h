#include <_.h>
#include "./bitmap.h"
#include <stdlib.h>
#include <algorithm>
#include <unordered_map>
#include <set>

struct simulate_result_16
{
    u_int16_t cycle, middle_product, nnz;
    u_int16_t mt125, mt25, mt50, mt75, mt100;
    u_int16_t nt125, nt25, nt50, nt75, nt100;
    u_int16_t C_ntr, A_ntr, B_ntr;
    double energy_cost, mac_rate, A_nut, B_nut, C_nut;

    simulate_result_16()
    {
        reset();
    }

    void reset(){
        cycle = middle_product = nnz = 0;
        A_nut = A_ntr = B_nut = B_ntr = C_nut = C_ntr = 0;
        nt125 = nt25 = nt50 = nt75 = nt100 = 0;
        mt125 = mt25 = mt50 = mt75 = mt100 = 0;
        energy_cost = 0;
    }

    // += 
    simulate_result_16& operator+=(const simulate_result_16 &a)
    {
        cycle += a.cycle;
        middle_product += a.middle_product;
        C_nut += a.C_nut;
        C_ntr += a.C_ntr;
        A_nut += a.A_nut;
        A_ntr += a.A_ntr;
        B_nut += a.B_nut;
        B_ntr += a.B_ntr;

        nt125 += a.nt125;
        nt25 += a.nt25;
        nt50 += a.nt50;
        nt75 += a.nt75;
        nt100 += a.nt100;
        
        mt125 += a.mt125;
        mt25 += a.mt25;
        mt50 += a.mt50;
        mt75 += a.mt75;
        mt100 += a.mt100;
        energy_cost += a.energy_cost;
        return *this;
    }
};

struct task_results
{
    u_int64_t ds_cycle, rm_cycle, smu_cycle, mp16, mp4, mpnz;

    u_int64_t ds_nt125, ds_nt25, ds_nt50, ds_nt75, ds_nt100;
    u_int64_t rm_nt125, rm_nt25, rm_nt50, rm_nt75, rm_nt100;
    u_int64_t smu_nt125, smu_nt25, smu_nt50, smu_nt75, smu_nt100;

    u_int64_t ds_mt125, ds_mt25, ds_mt50, ds_mt75, ds_mt100;
    u_int64_t rm_mt125, rm_mt25, rm_mt50, rm_mt75, rm_mt100;
    u_int64_t smu_mt125, smu_mt25, smu_mt50, smu_mt75, smu_mt100;

    u_int64_t ds_a_ntr, ds_b_ntr, ds_c_ntr;
    u_int64_t rm_a_ntr, rm_b_ntr, rm_c_ntr;
    u_int64_t smu_a_ntr, smu_b_ntr, smu_c_ntr;

    double ds_a_nut, ds_b_nut, ds_c_nut;
    double rm_a_nut, rm_b_nut, rm_c_nut;
    double smu_a_nut, smu_b_nut, smu_c_nut;

    double ds_mac, rm_mac, smu_mac;
    double ds_energy, rm_energy, smu_energy;
    double smu_speedup;

    task_results()
    {
        ds_cycle = rm_cycle = smu_cycle = 0;
        mp16 = mp4 = mpnz = 0;
        ds_mac = rm_mac = smu_mac = ds_energy = rm_energy = smu_energy = smu_speedup = 0;

        ds_nt125 = ds_nt25 = ds_nt50 = ds_nt75 = ds_nt100 = 0;
        rm_nt125 = rm_nt25 = rm_nt50 = rm_nt75 = rm_nt100 = 0;
        smu_nt125 = smu_nt25 = smu_nt50 = smu_nt75 = smu_nt100 = 0;

        ds_mt125 = ds_mt25 = ds_mt50 = ds_mt75 = ds_mt100 = 0;
        rm_mt125 = rm_mt25 = rm_mt50 = rm_mt75 = rm_mt100 = 0;
        smu_mt125 = smu_mt25 = smu_mt50 = smu_mt75 = smu_mt100 = 0;
        ds_a_ntr = ds_a_nut = ds_b_ntr = ds_b_nut = ds_c_ntr = ds_c_nut = 0;
        rm_a_ntr = rm_a_nut = rm_b_ntr = rm_b_nut = rm_c_ntr = rm_c_nut = 0;
        smu_a_ntr = smu_a_nut = smu_b_ntr = smu_b_nut = smu_c_ntr = smu_c_nut = 0;
    }

    task_results& operator+=(const task_results &a)
    {
        ds_cycle += a.ds_cycle;
        rm_cycle += a.rm_cycle;
        smu_cycle += a.smu_cycle;
        mp16 += a.mp16;
        mp4 += a.mp4;
        mpnz += a.mpnz;
        ds_energy += a.ds_energy;
        rm_energy += a.rm_energy;
        smu_energy += a.smu_energy;
        ds_a_ntr += a.ds_a_ntr, ds_a_nut += a.ds_a_nut, ds_b_ntr += a.ds_b_ntr, ds_b_nut += a.ds_b_nut, ds_c_ntr += a.ds_c_ntr, ds_c_nut += a.ds_c_nut;
        rm_a_ntr += a.rm_a_ntr, rm_a_nut += a.rm_a_nut, rm_b_ntr += a.rm_b_ntr, rm_b_nut += a.rm_b_nut, rm_c_ntr += a.rm_c_ntr, rm_c_nut += a.rm_c_nut;
        smu_a_ntr += a.smu_a_ntr, smu_a_nut += a.smu_a_nut, smu_b_ntr += a.smu_b_ntr, smu_b_nut += a.smu_b_nut, smu_c_ntr += a.smu_c_ntr, smu_c_nut += a.smu_c_nut;
        ds_nt125 += a.ds_nt125, ds_nt25 += a.ds_nt25, ds_nt50 += a.ds_nt50, ds_nt75 += a.ds_nt75, ds_nt100 += a.ds_nt100;
        rm_nt125 += a.rm_nt125, rm_nt25 += a.rm_nt25, rm_nt50 += a.rm_nt50, rm_nt75 += a.rm_nt75, rm_nt100 += a.rm_nt100;
        smu_nt125 += a.smu_nt125, smu_nt25 += a.smu_nt25, smu_nt50 += a.smu_nt50, smu_nt75 += a.smu_nt75, smu_nt100 += a.smu_nt100;
        ds_mt125 += a.ds_mt125, ds_mt25 += a.ds_mt25, ds_mt50 += a.ds_mt50, ds_mt75 += a.ds_mt75, ds_mt100 += a.ds_mt100;
        rm_mt125 += a.rm_mt125, rm_mt25 += a.rm_mt25, rm_mt50 += a.rm_mt50, rm_mt75 += a.rm_mt75, rm_mt100 += a.rm_mt100;
        smu_mt125 += a.smu_mt125, smu_mt25 += a.smu_mt25, smu_mt50 += a.smu_mt50, smu_mt75 += a.smu_mt75, smu_mt100 += a.smu_mt100;
        return *this;
    }

    void add_simulate_result(const simulate_result_16 &a, int stc_id, int times = 1)
    {
        switch (stc_id)
        {
        case 1: // DS-STC
            ds_cycle += a.cycle * times;
            ds_energy += a.energy_cost * times;
            
            ds_nt125 += a.nt125 * times, ds_nt25 += a.nt25 * times, ds_nt50 += a.nt50 * times, ds_nt75 += a.nt75 * times, ds_nt100 += a.nt100 * times;
            ds_mt125 += a.mt125 * times, ds_mt25 += a.mt25 * times, ds_mt50 += a.mt50 * times, ds_mt75 += a.mt75 * times, ds_mt100 += a.mt100 * times;
            ds_a_ntr += a.A_ntr * times, ds_a_nut += a.A_nut * times, ds_b_ntr += a.B_ntr * times, ds_b_nut += a.B_nut * times, ds_c_ntr += a.C_ntr * times, ds_c_nut += a.C_nut * times;
            break;
        case 2: // RM-STC
            rm_cycle += a.cycle * times;
            rm_energy += a.energy_cost * times;
            rm_nt125 += a.nt125 * times, rm_nt25 += a.nt25 * times, rm_nt50 += a.nt50 * times, rm_nt75 += a.nt75 * times, rm_nt100 += a.nt100 * times;
            rm_mt125 += a.mt125 * times, rm_mt25 += a.mt25 * times, rm_mt50 += a.mt50 * times, rm_mt75 += a.mt75 * times, rm_mt100 += a.mt100 * times;
            rm_a_ntr += a.A_ntr * times, rm_a_nut += a.A_nut * times, rm_b_ntr += a.B_ntr * times, rm_b_nut += a.B_nut * times, rm_c_ntr += a.C_ntr * times, rm_c_nut += a.C_nut * times;
            break;
        case 3: // SMU
            smu_cycle += a.cycle * times;
            smu_energy += a.energy_cost * times;
            smu_nt125 += a.nt125 * times, smu_nt25 += a.nt25 * times, smu_nt50 += a.nt50 * times, smu_nt75 += a.nt75 * times, smu_nt100 += a.nt100 * times;
            smu_mt125 += a.mt125 * times, smu_mt25 += a.mt25 * times, smu_mt50 += a.mt50 * times, smu_mt75 += a.mt75 * times, smu_mt100 += a.mt100 * times;
            smu_a_ntr += a.A_ntr * times, smu_a_nut += a.A_nut * times, smu_b_ntr += a.B_ntr * times, smu_b_nut += a.B_nut * times, smu_c_ntr += a.C_ntr * times, smu_c_nut += a.C_nut * times;
            mpnz += a.middle_product * times;
            break;
        default:
            break;
        }
    }
};

#pragma omp declare reduction(+ : task_results : omp_out += omp_in)\
    initializer(omp_priv = task_results())

struct smu_task
{
    TileIndex task_id;
    TileIndex task_mask;

    bool operator<(const smu_task &a) const
    {
        return TileBitmap_popcount(task_mask) > TileBitmap_popcount(a.task_mask);
    }
};

struct smu_scheduler_task
{
    TileIndex row, col, k, nnz;
    TileBitmap a, b, c;

    smu_scheduler_task()
    {
        row = col = k = nnz = 0;
        a = b = c = 0;
    }
};

struct smu_scheduler
{
    TileIndex task_num, total_nnz, busy;
    smu_scheduler_task task[16];

    smu_scheduler()
    {
        task_num = total_nnz = busy = 0;
    }
};

enum energy_model
{
    MAC_FP64, SMU_CTRL, SMU_SCHEDULER,
    LINE_BUFFER_WRITE, LINE_BUFFER_READ, // FP64  (0.5 KB)
    HALF_BUFFER_READ, HALF_BUFFER_WRITE, // FP64 （1KB）
    BUFFER_READ,   BUFFER_WRITE,         // FP64 （2KB）
    REG_U8_READ,   REG_U8_WRITE,         // U8
    REG_U16_READ,  REG_U16_WRITE,        // U16
    REG_FP64_READ, REG_FP64_WRITE,       // FP64  (64 KB)
    SMU_U8_READ,   SMU_U8_WRITE,         // U8
    SMU_U12_READ,  SMU_U12_WRITE,        // U12
    SMU_U16_READ,  SMU_U16_WRITE,        // U16
    DS_SCATTER,    DS_GATHER,            // DS
    RM_SCATTER,    RM_GATHER,            // RM
};

double energy_costs[] = {
    8.8014,  0.1368, 2.0944,
    0.69605, 0.69605,
    1.39211, 1.39211,
    2.59238, 2.78587,
    0.47821, 0.47821,
    0.95642, 0.95642,
    4.69205, 4.74013,
    0.15941, 0.15941,
    0.17024, 0.17024,
    0.18108, 0.18108,
    4.24768, 0.32836,
    1.06192, 1.06192
};

u_int16_t Tile_16_Bitmap_x_Bitmap(
    TileBitmap a16, TileBitmap *a4, // input A
    TileBitmap b16, TileBitmap *b4  // input B
)
{
    TileBitmap idx[4] = {0, 1, 17, 273}, nnzA = TileBitmap_popcount(a16);
    TileBitmap C_bitmap_ptr[16] = {0}, nnzCnt = 0;

    for (TileIndex Ai = 0; Ai < nnzA; ++Ai)
    {
        TileBitmap Ak = TileBitmap_lowbit(a16);
        a16 ^= Ak;

        TileIndex Actz = TileBitmap_ctz(Ak);
        TileBitmap cur_A_4_bitmap = a4[Ai];
        if (!cur_A_4_bitmap)
            continue;

        TileIndex row = Actz & 0xc, Acol = Actz & 3, Bj = TileBitmap_popcount(b16 & (0xf * idx[Acol]));
        TileIndex B_row_map = (b16 & (0xf << (Acol << 2))) >> (Acol << 2); // only using 4 bits

        while (B_row_map)
        {
            TileIndex Bk = TileBitmap_lowbit(B_row_map);
            B_row_map ^= Bk;

            TileIndex col = TileBitmap_ctz(Bk);

            TileBitmap cur_B_4_bitmap = b4[Bj++];
            if (!cur_B_4_bitmap)
                continue;

            C_bitmap_ptr[row | col] |= TileBitmap_x(cur_A_4_bitmap, cur_B_4_bitmap);
        }
    }

#pragma unroll
    for (TileIndex i = 0; i < 16; ++i)
    {
        nnzCnt += TileBitmap_popcount(C_bitmap_ptr[i]);
    }
    return nnzCnt;
}

template <int mac_num, int scheduler_num>
void SMX_SpMV_Task_Gen(smu_scheduler *s, TileBitmap a160, TileBitmap *a40, TileBitmap a161, TileBitmap *a41, double *control_energy, bool read_y, TileBitmap bv)
{
    TileIndex scheduler_num_mask = scheduler_num - 1, scheduler_num_bit = TileBitmap_ctz(scheduler_num);
    double energy = 0;
    smu_task mask[8] = {0};
    TileIndex cnt[8] = {0}, idx = 0, idx0 = 0;
    TileBitmap res = 0, _a160 = a160, _a161 = a161;
    TileIndex nnzA0 = TileBitmap_popcount(a160), nnzA1 = TileBitmap_popcount(a161);
    if (read_y) energy += 32 * energy_costs[HALF_BUFFER_READ];
    energy += 2 * energy_costs[SMU_CTRL];

    TileBitmap bv16 = 0, bv4[4] = {0};
    for (int i = 0; i < 4; ++i) {
        TileBitmap cur_bv = (bv >> (i << 2)) & 0xf;
        if (cur_bv) 
        {
            bv16 |= 1 << (i << 2);
            for (int j = 0; j < 4; ++j) if (cur_bv & (1 << j)) bv4[i] |= 1 << (j << 2);
        }
    }

#pragma unroll
    for (TileIndex i = 0; i < 4; ++i)
    {
        res = TileBitmap_col_x_realrow(a160 & 0x1111, bv4[i]);

        while (res)
        {
            TileBitmap k = res & -res;
            cnt[TileBitmap_ctz(k) >> 2] |= 1 << i;
            res ^= k;
        }

        a160 >>= 1;
    }

    if (a161)
    {
#pragma unroll
        for (TileIndex i = 0; i < 4; ++i)
        {
            res = TileBitmap_col_x_realrow(a161 & 0x1111, bv4[i]);

            while (res)
            {
                TileBitmap k = res & -res;
                cnt[TileBitmap_ctz(k) >> 2 | 4] |= 1 << i;
                res ^= k;
            }

            a161 >>= 1;
        }
    }

#pragma unroll
    for (TileIndex i = 0; i < 4; ++i)
    {
        if (cnt[i] == 0)
            continue;
        mask[idx].task_id = i;
        mask[idx].task_mask = cnt[i];
        ++idx;
    }
    idx0 = idx;

#pragma unroll
    for (TileIndex i = 4; i < 8; ++i)
    {
        if (cnt[i] == 0)
            continue;
        mask[idx].task_id = i;
        mask[idx].task_mask = cnt[i];
        ++idx;
    }

    for (TileIndex i = 0; i < idx0; ++i)
    {
        TileIndex row = mask[i].task_id;
        TileIndex _m = mask[i].task_mask;

        while (_m)
        {
            TileIndex k = _m & -_m, ctz = TileBitmap_ctz(k);
            _m ^= k;
            TileBitmap idxa = (1 << (row << 2 | ctz)) - 1;
            idxa = TileBitmap_popcount(idxa & _a160);
            // TileIndex _nnz = TileBitmap_popcount(a40[idxa]);
            TileIndex _nnz = TileBitmap_x_mpd(a40[idxa], bv4[ctz]);
            if (_nnz) {
                TileIndex _idx = s[i].task_num++;
                s[i].total_nnz += _nnz;
                s[i].task[_idx].row = row;
                s[i].task[_idx].col = ctz;
                s[i].task[_idx].k = ctz;
                s[i].task[_idx].nnz = _nnz;
                s[i].task[_idx].a = a40[idxa];
                s[i].task[_idx].b = bv4[ctz];
                s[i].task[_idx].c = TileBitmap_x(a40[idxa], bv4[ctz]);
                energy += energy_costs[SMU_U8_WRITE];
            }
        }
    }

    for (TileIndex i = idx0; i < idx; ++i)
    {
        TileIndex row = mask[i].task_id - 4;
        TileIndex _m = mask[i].task_mask;
        while (_m)
        {
            TileIndex k = _m & -_m, ctz = TileBitmap_ctz(k);
            _m ^= k;
            TileBitmap idxa = (1 << (row << 2 | ctz)) - 1;
            idxa = TileBitmap_popcount(idxa & _a161);
            // TileIndex _nnz = TileBitmap_popcount(a41[idxa]);
            TileIndex _nnz = TileBitmap_x_mpd(a41[idxa], bv4[ctz]);
            if (_nnz) {
                TileIndex id = (i - idx0) & scheduler_num_mask;
                TileIndex _idx = s[id].task_num++;
                s[id].total_nnz += _nnz;
                s[id].task[_idx].row = row;
                s[id].task[_idx].col = ctz;
                s[id].task[_idx].k = ctz;
                s[id].task[_idx].nnz = _nnz;
                s[id].task[_idx].a = a41[idxa];
                s[id].task[_idx].b = bv4[ctz];
                s[id].task[_idx].c = TileBitmap_x(a41[idxa], bv4[ctz]);
                energy += energy_costs[SMU_U8_WRITE];
            }
        }
    }

    *control_energy = energy;
}

u_int16_t nnzCount(TileBitmap x16, TileBitmap *x4)
{
    u_int16_t nnz = 0;
    TileIndex nnzX = TileBitmap_popcount(x16);
    for (TileIndex i = 0; i < nnzX; ++i)
    {
        nnz += TileBitmap_popcount(x4[i]);
    }
    return nnz;
}

template <int mac_num, int scheduler_num>
void SMX_SpGEMM_Task_Gen(smu_scheduler *s, TileBitmap a16, TileBitmap *a4, TileBitmap b16, TileBitmap *b4, double *control_energy, bool load_a)
{
    TileIndex scheduler_num_mask = scheduler_num - 1, scheduler_num_bit = TileBitmap_ctz(scheduler_num);
    double energy = 0;
    smu_task mask[16] = {0};
    TileIndex cnt[16] = {0}, idx = 0;
    TileBitmap res = 0, _a16 = a16, _b16 = b16;
    u_int16_t nnzAnum = nnzCount(a16, a4), nnzCnum = Tile_16_Bitmap_x_Bitmap(a16, a4, b16, b4);
    energy_model buffer_itype = scheduler_num == 8? HALF_BUFFER_READ : BUFFER_READ;
    energy_model buffer_otype = scheduler_num == 8? HALF_BUFFER_WRITE : BUFFER_WRITE;

    if (load_a) energy += nnzAnum * (energy_costs[BUFFER_WRITE] + energy_costs[REG_FP64_READ]);
    energy += nnzCnum * (energy_costs[buffer_otype] + energy_costs[buffer_itype]);
    energy += energy_costs[SMU_CTRL];
    
    int8_t sid = 0, forward = 1;
    
#pragma unroll
    for (TileIndex i = 0; i < 4; ++i)
    {
        TileBitmap cur_a = a16 & 0x1111, cur_b = b16 & 0xf;
        
        res = TileBitmap_col_x_row(cur_a, cur_b);
        cur_a = TileBitmap_popcount(cur_a);
        cur_b = TileBitmap_popcount(cur_b);

        if (cur_a <= cur_b) {
            for (int k = 0; k < 4; ++k) for (int j = 0; j < 4; ++j) 
            {
                if (res & (1 << (j << 2 | k)))
                {
                    TileIndex row = j, col = k;
                    TileBitmap idxa = (1 << (row << 2 | i)) - 1,
                        idxb = (1 << (i << 2 | col)) - 1;

                    idxa = TileBitmap_popcount(idxa & _a16);
                    idxb = TileBitmap_popcount(idxb & _b16);
                    TileIndex _nnz = TileBitmap_x_mpd(a4[idxa], b4[idxb]);
                    TileIndex _idx = s[sid].task_num++;
                    s[sid].total_nnz += _nnz;
                    s[sid].task[_idx].row = row;
                    s[sid].task[_idx].col = col;
                    s[sid].task[_idx].k = i;
                    s[sid].task[_idx].nnz = _nnz;
                    s[sid].task[_idx].a = a4[idxa];
                    s[sid].task[_idx].b = b4[idxb];
                    s[sid].task[_idx].c = TileBitmap_x(a4[idxa], b4[idxb]);
                    energy += energy_costs[SMU_U8_WRITE];
                    sid = (sid + 1) & scheduler_num_mask;
                }
            }
        } else {
            while (res)
            {
                TileBitmap k = res & -res;
                TileIndex ctz = TileBitmap_ctz(k), row = ctz >> 2, col = ctz & 3;
                TileBitmap idxa = (1 << (row << 2 | i)) - 1,
                        idxb = (1 << (i << 2 | col)) - 1;

                idxa = TileBitmap_popcount(idxa & _a16);
                idxb = TileBitmap_popcount(idxb & _b16);
                TileIndex _nnz = TileBitmap_x_mpd(a4[idxa], b4[idxb]);
                TileIndex _idx = s[sid].task_num++;
                s[sid].total_nnz += _nnz;
                s[sid].task[_idx].row = row;
                s[sid].task[_idx].col = col;
                s[sid].task[_idx].k = i;
                s[sid].task[_idx].nnz = _nnz;
                s[sid].task[_idx].a = a4[idxa];
                s[sid].task[_idx].b = b4[idxb];
                s[sid].task[_idx].c = TileBitmap_x(a4[idxa], b4[idxb]);
                energy += energy_costs[SMU_U8_WRITE];
                sid = (sid + 1) & scheduler_num_mask;
                // sid += forward;
                // if (forward == 1 && sid == scheduler_num - 1) forward = -1;
                // else if (forward == -1 && sid == 0) forward = 1;
                res ^= k;
            }
        }

        a16 >>= 1;
        b16 >>= 4;
    }

    *control_energy = energy;
}

template <int mac_num, int scheduler_num>
simulate_result_16 SMU_scheduler_simulate(smu_scheduler *s)
{
    TileIndex scheduler_num_mask = scheduler_num - 1;
    TileIndex cursor[scheduler_num] = {0}, cycle = 0, ci = 0;
    smu_scheduler_task dot4buffer[scheduler_num];
    u_int16_t cur_mac = 0;
    TileIndex flag = 1;
    double energy = 0;
    simulate_result_16 res;

    for (int si = 0; si < scheduler_num; ++si) {
        for (int ti = 0; ti < s[si].task_num; ++ti) {
            u_int8_t cnz = TileBitmap_popcount(s[si].task[ti].c);
            res.C_ntr += cnz;
            double network_util = cnz * 1.0 / 16 * 100;
            if (network_util <= 12.5) res.nt125 += cnz;
            else if (network_util <= 25) res.nt25 += cnz;
            else if (network_util <= 50) res.nt50 += cnz;
            else if (network_util <= 75) res.nt75 += cnz;
            else res.nt100 += cnz;
        }
    }

    while (1)
    {
        flag = 1;
        if (cycle)
        {
            std::unordered_map<TileIndex, TileIndex> A_blocks, B_blocks;
            TileIndex consum = 0;

            TileIndex it = 0, executed = 0, c_ntr = 0, a_ntr = 0, b_ntr = 0;
            ci = cycle & 1? 0: scheduler_num - 1;
            TileBitmap occupy = 0, cb = 0;
            while (it < scheduler_num && consum < mac_num)
            {
                if (s[ci].busy)
                {
                    cb = 1 << (dot4buffer[ci].row << 2 | dot4buffer[ci].col);
                    if ((consum < mac_num || dot4buffer[ci].nnz == 0) && !(occupy & cb))
                    {
                        occupy |= cb;
                        if (consum + dot4buffer[ci].nnz <= mac_num)
                        {
                            s[ci].busy = 0;
                            consum += dot4buffer[ci].nnz;
                            A_blocks[dot4buffer[ci].row << 2 | dot4buffer[ci].k] = TileBitmap_popcount(dot4buffer[ci].a);
                            #ifdef SPMV
                                energy += 2.8 * (energy_costs[HALF_BUFFER_READ] + energy_costs[HALF_BUFFER_WRITE]);
                                B_blocks[dot4buffer[ci].k] = 4;
                            #else
                                energy += dot4buffer[ci].nnz * (energy_costs[LINE_BUFFER_READ] + energy_costs[LINE_BUFFER_WRITE]);
                                B_blocks[dot4buffer[ci].k << 2 | dot4buffer[ci].col] = TileBitmap_popcount(dot4buffer[ci].b);
                            #endif
                            ++executed;
                            c_ntr += TileBitmap_popcount(dot4buffer[ci].c);
                        }
                        else
                        {
                            TileIndex deal_nnz = mac_num - consum;
                            dot4buffer[ci].nnz -= deal_nnz;
                            consum = mac_num;
                            #ifdef SPMV
                                energy += 2.8 * (energy_costs[HALF_BUFFER_READ] + energy_costs[HALF_BUFFER_WRITE]);
                            #else
                                energy += deal_nnz * (energy_costs[LINE_BUFFER_READ] + energy_costs[LINE_BUFFER_WRITE]);
                            #endif
                            // A_blocks[dot4buffer[ci].row << 2 | dot4buffer[ci].k] = TileBitmap_popcount(dot4buffer[ci].a);
                            // B_blocks[dot4buffer[ci].k << 2 | dot4buffer[ci].col] = TileBitmap_popcount(dot4buffer[ci].b);
                            ++executed;
                            c_ntr += (deal_nnz + 3) / 4;
                        }
                    }
                }
                ci = cycle & 1? ci + 1: ci - 1;
                ++it;
            }
            res.C_nut += sqrt(executed * 1.0 / scheduler_num) * c_ntr;
            
            for (auto &a : A_blocks)
            {
                energy += a.second * energy_costs[BUFFER_READ];
                a_ntr += a.second;
            }
            for (auto &b : B_blocks)
            {
                energy += b.second * energy_costs[BUFFER_READ];
                b_ntr += b.second;
            }
            res.A_ntr += a_ntr;
            res.B_ntr += b_ntr;
            res.A_nut += sqrt(executed / 16.0 / scheduler_num) * a_ntr;
            res.B_nut += sqrt(executed / 16.0 / scheduler_num) * b_ntr;
            double mac_rate = consum * 100.0 / mac_num;
            if (mac_rate <= 12.5) res.mt125 += 1;
            else if (mac_rate <= 25) res.mt25 += 1;
            else if (mac_rate <= 50) res.mt50 += 1;
            else if (mac_rate <= 75) res.mt75 += 1;
            else res.mt100 += 1;
        }

        bool enq_flag = false;

        for (TileIndex si = 0; si < scheduler_num; ++si)
        {
            if (cursor[si] < s[si].task_num)
            {
                if (!s[si].busy)
                {
                    while (!s[si].task[cursor[si]].nnz) {
                        ++cursor[si];
                        if (cursor[si] >= s[si].task_num) break;
                    }
                    if (s[si].task[cursor[si]].nnz) {
                        cur_mac += s[si].task[cursor[si]].nnz;
                        dot4buffer[si] = s[si].task[cursor[si]];
                        energy += energy_costs[SMU_U8_READ] + s[si].task[cursor[si]].nnz * energy_costs[SMU_U12_WRITE];
                        ++cursor[si];
                        s[si].busy = 1;
                        enq_flag = true;
                    }
                }
            }
            if (s[si].busy)
                flag = 0;
        }
        if (enq_flag) energy += energy_costs[SMU_SCHEDULER] * (scheduler_num >> 3);
        ++cycle;
        if (flag)
        {
            break;
        }
    }

    if (cycle) cycle -= 1;
    res.cycle = cycle;
    res.energy_cost = energy;
    res.middle_product = cur_mac;
    return res;
}

template <int mac_num, int scheduler_num>
simulate_result_16 SMU_SpMV(TileBitmap a160, TileBitmap *a40, TileBitmap a161, TileBitmap *a41, bool read_y, TileBitmap bv)
{
    double control_energy = 0, mac_energy = 0;
    smu_scheduler s[scheduler_num];
    SMX_SpMV_Task_Gen<mac_num, scheduler_num>(s, a160, a40, a161, a41, &control_energy, read_y, bv);
    simulate_result_16 res = SMU_scheduler_simulate<mac_num, scheduler_num>(s);
    res.energy_cost = (control_energy + res.energy_cost) / 1000; // (mJ)
    return res;
}

template <int mac_num, int scheduler_num>
simulate_result_16 SMU_SpGEMM(TileBitmap a16, TileBitmap *a4, TileBitmap b16, TileBitmap *b4, bool load_a)
{
    double control_energy = 0, mac_energy = 0;
    smu_scheduler s[scheduler_num];
    SMX_SpGEMM_Task_Gen<mac_num, scheduler_num>(s, a16, a4, b16, b4, &control_energy, load_a);
    simulate_result_16 res = SMU_scheduler_simulate<mac_num, scheduler_num>(s);
    res.energy_cost = (control_energy + res.energy_cost) / 1000; // (mJ)
    return res;
}

template <int mac_num>
simulate_result_16 Tile_16_Mx_ds(
    TileBitmap a16, TileBitmap *a4,
    TileBitmap b16, TileBitmap *b4
)
{
    simulate_result_16 res;
    double energy = 0;
    u_int8_t cycle = 0;
    TileBitmap amask[4] = {0x1111, 0x2222, 0x4444, 0x8888};
    TileBitmap bmask[4] = {0xf, 0xf0, 0xf00, 0xf000};
    TileIndex mac_col_len = mac_num >> 3;

    energy += Tile_16_Bitmap_x_Bitmap(a16, a4, b16, b4) * (energy_costs[BUFFER_WRITE] + energy_costs[BUFFER_READ]);

    for (TileIndex i = 0; i < 4; ++i)
    {
        TileBitmap A_col = a16 & amask[i], B_row = b16 & bmask[i];
        if (A_col && B_row)
        {
            for (TileIndex j = 0; j < 4; ++j)
            {
                TileBitmap _a = A_col, _b = B_row;
                u_int16_t A_col_len = 0;
                u_int16_t B_row_len = 0;

                while (_a)
                {
                    TileBitmap k = TileBitmap_lowbit(_a);
                    _a ^= k;

                    TileIndex ctz = TileBitmap_ctz(k);
                    TileBitmap a4b = a4[TileBitmap_popcount(a16 & ((1 << ctz) - 1))];
                    A_col_len += TileBitmap_popcount(a4b & amask[j]);
                }

                while (_b)
                {
                    TileBitmap k = TileBitmap_lowbit(_b);
                    _b ^= k;

                    TileIndex ctz = TileBitmap_ctz(k);
                    TileBitmap b4b = b4[TileBitmap_popcount(b16 & ((1 << ctz) - 1))];
                    B_row_len += TileBitmap_popcount(b4b & bmask[j]);
                }

                if (A_col_len && B_row_len)
                {
                    energy += (A_col_len + B_row_len) * (energy_costs[LINE_BUFFER_WRITE] + energy_costs[LINE_BUFFER_READ] + energy_costs[REG_FP64_READ]);
                    energy += A_col_len * B_row_len * (energy_costs[BUFFER_WRITE] + energy_costs[BUFFER_READ]);
                    energy += energy_costs[DS_GATHER] + energy_costs[DS_SCATTER]; // 数据聚集
                    for (int ai = 0; ai < A_col_len; ai += 8) {
                        for (int bi = 0; bi < B_row_len; bi += mac_col_len){
                            int cur_A_len = std::min(8, A_col_len - ai);
                            int cur_B_len = std::min((int)mac_col_len, B_row_len - bi);
                            int intermedia = cur_A_len * cur_B_len;
                            double mac_rate = intermedia * 100.0 / mac_num;
                            double network_rate = intermedia * 1.0 / 256 * 100;
                            if (mac_rate <= 12.5) res.mt125 += 1;
                            else if (mac_rate <= 25) res.mt25 += 1;
                            else if (mac_rate <= 50) res.mt50 += 1;
                            else if (mac_rate <= 75) res.mt75 += 1;
                            else res.mt100 += 1;
                            if (network_rate <= 12.5) res.nt125 += intermedia;
                            else if (network_rate <= 25) res.nt25 += intermedia;
                            else if (network_rate <= 50) res.nt50 += intermedia;
                            else if (network_rate <= 75) res.nt75 += intermedia;
                            else res.nt100 += intermedia;
                        }
                    }
                    res.A_ntr += A_col_len * ((B_row_len + mac_col_len - 1) / mac_col_len);
                    res.B_ntr += B_row_len * ((A_col_len + 7) >> 3);
                    res.C_ntr += A_col_len * B_row_len;
                    A_col_len = (A_col_len + 7) >> 3;
                    B_row_len = (B_row_len + mac_col_len - 1) / mac_col_len;
                    cycle += std::max(A_col_len * B_row_len, 1);
                }
            }
        }
    }

    res.cycle = res.A_nut = res.B_nut = res.C_nut = cycle;
    res.energy_cost = energy / 1000;
    return res;
}

template <int mac_num>
simulate_result_16 Tile_16_Mx_rm(
    TileBitmap a16, TileBitmap *a4,
    TileBitmap b16, TileBitmap *b4
)
{
    simulate_result_16 res;
    double energy = 0;
    u_int8_t cycle = 0;
    u_int16_t Arows[17] = {0}, Brows[17] = {0}, nzA = 0, nzB = 0;
    u_int8_t Acols[256] = {0}, Bcols[256] = {0};
    TileIndex mac_col_len, i_max, j_max, i_range;
    if (mac_num == 128) {
        mac_col_len = 4; // [4, 8, 16]
        j_max = 128 / mac_col_len / 2;
        i_max = 16 / j_max;
        i_range = 16;
    } else {
        mac_col_len = 4; // [4, 8]
        j_max = 64 / mac_col_len / 2;
        i_max = 16 / j_max;
        i_range = 8;
    }

    energy += Tile_16_Bitmap_x_Bitmap(a16, a4, b16, b4) * (energy_costs[LINE_BUFFER_WRITE] + energy_costs[LINE_BUFFER_READ]);
    energy += nnzCount(a16, a4) * (energy_costs[REG_FP64_READ]);

    bitmap2csr(a16, a4, Arows, Acols, &nzA);
    bitmap2csr(b16, b4, Brows, Bcols, &nzB);

    // std::set<u_int8_t> j_set;
    std::set<u_int8_t> col_set;

    for (u_int8_t i = 0; i < i_max; ++i)
    {
        u_int8_t rowlen[16], ib = i * j_max;

        for (u_int8_t j = 0; j < j_max; ++j)
            rowlen[j] = Arows[ib + j + 1] - Arows[ib + j];

        u_int8_t max_row_len = 0; 
        for (u_int8_t j = 0; j < j_max; ++j)
            max_row_len = std::max(max_row_len, rowlen[j]);

        for (u_int8_t window_base = 0; window_base < max_row_len; window_base += 2)
        {
            std::unordered_map<u_int8_t, u_int8_t> col_len, col_cnt;
            // j_set.clear();
            u_int16_t max_col_len = 0, c_ntr[4] = {0}, c_nut[4] = {0};
            
            for (u_int8_t j = 0; j < j_max; ++j)
            {
                u_int16_t Brow[2][16], max_len = 0;
                memset(Brow[0], -1, sizeof(Brow[0]));
                memset(Brow[1], -1, sizeof(Brow[1]));
                for (u_int8_t si = 0; si < 2; ++si) {
                    if (window_base + si < rowlen[j])
                    {
                        u_int8_t col = Acols[Arows[i * i_range + j] + window_base + si];
                        u_int16_t len = Brows[col + 1] - Brows[col];
                        max_col_len = std::max(max_col_len, len);
                        energy += (energy_costs[LINE_BUFFER_WRITE] + energy_costs[LINE_BUFFER_READ]) * len;// + len * energy_costs[MAC_FP64];
                        col_len[col] = len;
                        col_cnt[col]++;
                        // j_set.insert(j);
                        for (int k = 0; k < len; ++k)
                        {
                            Brow[si][k] = Bcols[Brows[col] + k];
                        }
                        max_len = std::max(max_len, len);
                    }
                }
                for (int ci = 0; ci < max_len; ci += mac_col_len)
                {
                    col_set.clear();
                    for (int si = 0; si < 2; ++si)
                    {
                        for (int cj = 0; cj < mac_col_len; ++cj)
                        {
                            if (Brow[si][ci + cj] != -1)
                            {
                                col_set.insert(Brow[si][ci + cj]);
                            }
                        }
                    }
                    u_int16_t cnz = col_set.size();
                    c_ntr[ci >> 2] += cnz;
                    c_nut[ci >> 2] += 1;
                    // c_ntr += cnz;
                    // c_nut += 1;
                    
                    double network_util = cnz * 1.0 / 16 * 100;
                    if (network_util <= 12.5) res.nt125 += cnz;
                    else if (network_util <= 25) res.nt25 += cnz;
                    else if (network_util <= 50) res.nt50 += cnz;
                    else if (network_util <= 75) res.nt75 += cnz;
                    else res.nt100 += cnz;
                }
            }
            res.C_ntr += c_ntr[0] + c_ntr[1] + c_ntr[2] + c_ntr[3];
            res.C_nut += sqrt(c_nut[0] * 1.0 / j_max) * c_ntr[0] + sqrt(c_nut[1] * 1.0 / j_max) * c_ntr[1] + sqrt(c_nut[2] * 1.0 / j_max) * c_ntr[2] + sqrt(c_nut[3] * 1.0 / j_max) * c_ntr[3];
            
            for (auto&it: col_len) {
                // res.B_ntr += it.second;
                energy += std::min(1.5, (double)col_cnt[it.first]) * it.second * (energy_costs[REG_FP64_READ] + energy_costs[HALF_BUFFER_WRITE]);
            }
            if (max_col_len)
                energy += energy_costs[RM_GATHER] + energy_costs[RM_SCATTER];

            int using_cycle = (max_col_len + mac_col_len - 1) / mac_col_len;
            // res.C_nut += j_set.size();
            cycle += using_cycle;

            for (int ci = 0; ci < using_cycle; ++ci)
            {
                int mac_sum = 0, len_start = ci * mac_col_len;
                for (auto&it: col_len) {
                    if (len_start < it.second) {
                        res.A_ntr += col_cnt[it.first];
                        res.B_ntr += std::min(it.second - len_start, (int)mac_col_len);
                        mac_sum += std::min((int)mac_col_len, it.second - len_start) * col_cnt[it.first];
                    }
                }
                // res.C_ntr += mac_sum;
                double mac_rate = mac_sum * 100.0 / (double)mac_num;
                if (mac_rate <= 12.5) res.mt125 += 1;
                else if (mac_rate <= 25) res.mt25 += 1;
                else if (mac_rate <= 50) res.mt50 += 1;
                else if (mac_rate <= 75) res.mt75 += 1;
                else res.mt100 += 1;
            }
        }
    }
    // res.A_ntr = nzA;
    res.cycle = res.A_nut = res.B_nut = cycle;
    res.energy_cost = energy / 1000;
    return res;
}

template <int mac>
simulate_result_16 Tile_16_Mx_smu(
    TileBitmap A_16_bitmap, TileBitmap *A_bitmap_ptr,
    TileBitmap B_16_bitmap, TileBitmap *B_bitmap_ptr, bool&load_a
)
{
    simulate_result_16 res = SMU_SpGEMM<mac, SCHEDULER>(A_16_bitmap, A_bitmap_ptr, B_16_bitmap, B_bitmap_ptr, load_a);
    load_a = false;
    return res;
}

template <int mac_num>
void Tile_16_M_x_M(
    TileBitmap A_16_bitmap, TileBitmap *A_bitmap_ptr,
    TileBitmap B_16_bitmap, TileBitmap *B_bitmap_ptr,
    // TileBitmap C_16_Bitmap, TileBitmap *C_bitmap_ptr,
    simulate_result_16 &_ds, simulate_result_16 &_rm, simulate_result_16 &_smu, bool&load_a)
{
    _smu = Tile_16_Mx_smu<mac_num>(A_16_bitmap, A_bitmap_ptr, B_16_bitmap, B_bitmap_ptr, load_a);
    _ds = Tile_16_Mx_ds<mac_num>(A_16_bitmap, A_bitmap_ptr, B_16_bitmap, B_bitmap_ptr);
    _rm = Tile_16_Mx_rm<mac_num>(A_16_bitmap, A_bitmap_ptr, B_16_bitmap, B_bitmap_ptr);
}

template <int mac_num>
simulate_result_16 TileMatrixSpMV_nv(TileBitmap A16b, TileBitmap *A4b)
{
    // return std::make_tuple(0, 4096 / mac_num);
    simulate_result_16 res;
    res.cycle = 4096 / mac_num;
    return res;
}

template <int mac_num>
simulate_result_16 TileMatrixSpMSpV_ds(TileBitmap A16b, TileBitmap *A4b, TileBitmap bv)
{
    simulate_result_16 res;
    double energy = 0;
    u_int8_t cycle = 0;
    TileBitmap amask[4] = {0x1111, 0x2222, 0x4444, 0x8888};
    TileIndex nnzA = TileBitmap_popcount(A16b);
    // energy += (1 + nnzA) * (energy_costs[REG_U16_READ] + energy_costs[SMU_U16_WRITE]);
    // energy += 16 * (energy_costs[REG_FP64_READ] + energy_costs[BUFFER_WRITE] + 2 * energy_costs[BUFFER_READ] + energy_costs[REG_FP64_WRITE]);
    energy += 16 * energy_costs[BUFFER_READ];

    for (TileIndex i = 0; i < 4; ++i)
    {
        TileBitmap A_col = A16b & amask[i];
        if (A_col)
        {
            for (TileIndex j = 0; j < 4; ++j)
            {
                TileBitmap _a = A_col;
                u_int16_t A_col_len = 0, B_row = bv & (1 << (i << 2 | j));
                if (!B_row) continue;

                while (_a)
                {
                    TileBitmap k = TileBitmap_lowbit(_a);
                    _a ^= k;

                    TileIndex ctz = TileBitmap_ctz(k);
                    TileBitmap a4b = A4b[TileBitmap_popcount(A16b & ((1 << ctz) - 1))];
                    A_col_len += TileBitmap_popcount(a4b & amask[j]);
                }

                if (A_col_len)
                {
                    energy += (A_col_len + 1) * (energy_costs[REG_FP64_READ] + energy_costs[LINE_BUFFER_WRITE] + energy_costs[LINE_BUFFER_READ]);
                    energy += A_col_len * (energy_costs[BUFFER_WRITE] + energy_costs[BUFFER_READ]);
                    energy += energy_costs[DS_GATHER] + energy_costs[DS_SCATTER]; // 数据聚集
                    for (int ai = 0; ai < A_col_len; ai += 8) {
                        int cur_A_len = std::min(8, A_col_len - ai);
                        int intermedia = cur_A_len;
                        double mac_rate = intermedia * 100.0 / mac_num;
                        double network_rate = intermedia * 1.0 / 256 * 100;
                        if (mac_rate <= 12.5) res.mt125 += 1;
                        else if (mac_rate <= 25) res.mt25 += 1;
                        else if (mac_rate <= 50) res.mt50 += 1;
                        else if (mac_rate <= 75) res.mt75 += 1;
                        else res.mt100 += 1;
                        if (network_rate <= 12.5) res.nt125 += intermedia;
                        else if (network_rate <= 25) res.nt25 += intermedia;
                        else if (network_rate <= 50) res.nt50 += intermedia;
                        else if (network_rate <= 75) res.nt75 += intermedia;
                        else res.nt100 += intermedia;
                    }
                    res.A_ntr += A_col_len;
                    res.B_ntr += ((A_col_len + 7) >> 3);
                    res.C_ntr += A_col_len;
                    A_col_len = (A_col_len + 7) >> 3;
                    cycle += std::max(A_col_len, (u_int16_t)1);
                }
            }
        }
    }
    // return std::make_tuple(midx, cycle);
    // return std::make_tuple(0, cycle);
    res.cycle = cycle;
    res.C_nut = cycle;
    res.energy_cost = energy / 1000;
    return res;
}

template <int mac_num>
simulate_result_16 TileMatrixSpMSpV_rm(TileBitmap A16b, TileBitmap *A4b, TileBitmap bv)
{
    simulate_result_16 res;
    double energy = 0;
    u_int16_t Arows[17] = {0}, nzA = 0;
    u_int8_t Acols[256] = {0}, cycle = 0, stride = 2;
    TileIndex mac_col_len, i_max, j_max;
    if (mac_num == 128) {
        mac_col_len = 4; // [4, 8, 16]
        j_max = 128 / mac_col_len / stride;
        i_max = 16 / j_max;
    } else {
        mac_col_len = 4; // [4, 8]
        j_max = 64 / mac_col_len / stride;
        i_max = 16 / j_max;
    }
    energy += 16 * energy_costs[LINE_BUFFER_READ];
    energy += nnzCount(A16b, A4b) * energy_costs[REG_FP64_READ];

    bitmap2csr(A16b, A4b, Arows, Acols, &nzA);
    for (u_int8_t i = 0; i < i_max; ++i)
    {
        u_int8_t rowlen[16], ib = i * j_max;
        for (u_int8_t j = 0; j < j_max; ++j)
            rowlen[j] = Arows[ib + j + 1] - Arows[ib + j];
        u_int8_t max_row_len = 0; 
        for (u_int8_t j = 0; j < j_max; ++j)
            max_row_len = std::max(max_row_len, rowlen[j]);
        for (u_int8_t window_base = 0; window_base < max_row_len; window_base += stride)
        {   
            std::unordered_map<u_int8_t, u_int8_t> col_len, col_cnt;
            int max_col_len = 0, A_elements = 0;
            double tmp_energy = 0;
            
            for (u_int8_t j = 0; j < j_max; ++j)
            {
                int c_ntr = 0;
                for (u_int8_t si = 0; si < stride; ++si)
                    if (window_base + si < rowlen[j])
                    {
                        u_int8_t col = Acols[Arows[i << 2 | j] + window_base + si];
                        if (bv & (1 << col)) {
                            tmp_energy += energy_costs[LINE_BUFFER_WRITE] + energy_costs[LINE_BUFFER_READ];
                            col_len[col] = 1;
                            col_cnt[col]++;
                            max_col_len = 1;
                            A_elements++;
                            c_ntr = 1;
                        }
                    }
                res.C_ntr += c_ntr;
                if (c_ntr) res.nt125 += 1;
            }
            int mac_sum = 0;

            for (auto&it: col_len) {
                tmp_energy += std::min(1.5, (double)col_cnt[it.first]) * it.second * (energy_costs[REG_FP64_READ] + energy_costs[HALF_BUFFER_WRITE]);
                res.A_ntr += col_cnt[it.first];
                res.B_ntr += 1;
                mac_sum += col_cnt[it.first];
            }

            if (max_col_len)
                tmp_energy += energy_costs[RM_GATHER] + energy_costs[RM_SCATTER];
            energy += tmp_energy;

            // res.C_ntr += mac_sum;
            double mac_rate = mac_sum * 100.0 / (double)mac_num;
            if (mac_rate <= 12.5) res.mt125 += 1;
            else if (mac_rate <= 25) res.mt25 += 1;
            else if (mac_rate <= 50) res.mt50 += 1;
            else if (mac_rate <= 75) res.mt75 += 1;
            else res.mt100 += 1;
        }
        cycle += (max_row_len + 1) >> 1;
    }
    res.cycle = cycle;
    res.energy_cost = energy / 1000;
    return res;
}

template <int mac_num>
simulate_result_16 TileMatrixSpMSpV_smu(TileBitmap A16b0, TileBitmap *A4b0, TileBitmap A16b1, TileBitmap *A4b1, bool read_y, TileBitmap bv)
{
    return SMU_SpMV<mac_num, 8>(A16b0, A4b0, A16b1, A4b1, read_y, bv);
}

template <int mac_num>
inline void TileMatrixSpMSpV_16(
    TileBitmap A16b0, TileBitmap *A4b0, TileBitmap A16b1, TileBitmap *A4b1,
    simulate_result_16 &_ds, simulate_result_16 &_rm, simulate_result_16 &_smu, bool read_y, TileBitmap bv = 0xffff)
{
    _ds = TileMatrixSpMSpV_ds<mac_num>(A16b0, A4b0, bv);
    _rm = TileMatrixSpMSpV_rm<mac_num>(A16b0, A4b0, bv);
    if (A16b1)
    {
        _ds += TileMatrixSpMSpV_ds<mac_num>(A16b1, A4b1, bv);
        _rm += TileMatrixSpMSpV_rm<mac_num>(A16b1, A4b1, bv);
    }
    _smu = TileMatrixSpMSpV_smu<mac_num>(A16b0, A4b0, A16b1, A4b1, read_y, bv);
}