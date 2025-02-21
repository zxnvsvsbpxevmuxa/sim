#pragma once
#include <tmatrix/Utils/timer.h>
#include <tmatrix/MMIO/mmio_highlevel.h>
#include <tmatrix/DataStructure/Factory.h>
#include <tmatrix/DataStructure/transform.h>
#include <tmatrix/Calculation/TileMatrix_x.h>

TileMatrix* load_matrix(const char*filepath)
{
    MatIndex m, n, nnz;
    MatIndex *row_ptr, *col_ptr;
    MatValue *val_ptr;
    Timer read_mtx, csr_to_seqtile;

    char tile_path[128];
    sprintf(tile_path, "%s.tile", filepath);
    bool exists = file_exists(tile_path);

    if (!exists) {
        int isSymmetric;
        echo(start_status, "Reading Matrix file: [bold green]%s[/]", filepath);
        timer_start(read_mtx);
        int status = mmio_allinone(filepath, &m, &n, &nnz, &isSymmetric, &row_ptr, &col_ptr, &val_ptr, true);
        timer_end(read_mtx);
        echo(stop_status, "");
        if (status)
        {
            echo(error, "Matrix is not square, please use square matrix");
            return nullptr;
        }
        echo(success, "m = %d, n = %d, nnz = %d; time: %f ms", m, n, nnz, timer_duration(read_mtx));
    }

    echo(start_status, "Converting CSR to TileMatrix");
    timer_start(csr_to_seqtile);
    TileMatrix *seqtile = CSR_to_TileMatrix(m, n, nnz, row_ptr, col_ptr, val_ptr, filepath, true, false, true);
    if (seqtile == nullptr)
    {
        echo(error, "Matrix is not symmetric, please use symmetric matrix");
        echo(stop_status, "");
        return nullptr;
    }

    symbolRateCalculation(seqtile);

    timer_end(csr_to_seqtile);
    echo(stop_status, "");
    echo(success, "Convert CSR to TileMatrix time: %f ms, m: %d, _nnz: %d, l1nnz: %d, nnz: %d", timer_duration(csr_to_seqtile), seqtile->meta_m, seqtile->_nnz, seqtile->tile_val_map[seqtile->_nnz], seqtile->meta_nnz);

    if (!exists) {
        free(row_ptr);
        free(col_ptr);
        free(val_ptr);
    }

    return seqtile;
}

int test_SpGEMM(const char *filepath, const char*filepath2 = nullptr)
{
    TileMatrix* seqtile_A = load_matrix(filepath), *seqtile_B;
    if (seqtile_A == nullptr) return 1;
    if (filepath2 != nullptr) {
        seqtile_B = load_matrix(filepath2);
        if (seqtile_B == nullptr) return 1;
    } else seqtile_B = seqtile_A;
    
    echo(start_status, "TMX Symbolic GEMM");
    TileMatrix *seqtile_c = TileMatrix_GEMM_Sym_Sim(seqtile_A, seqtile_B);
    // echo(info, "TMX Sym: _nnz: %d, l1nnz: %d, nnz: %d", seqtile_c->_nnz, seqtile_c->tile_val_map[seqtile_c->_nnz], seqtile_c->meta_nnz);
    echo(info, "TMX Sym: _nnz: %d", seqtile_c->_nnz);
    // symbolRateCalculation(seqtile_c);
    echo(start_status, "Running Simulation 128 MAC...");
    task_results t32 = TileMatrix_GEMM_Num_Sim<128>(seqtile_A, seqtile_B, seqtile_c), t64;
    if (seqtile_B == seqtile_A) {
        echo(start_status, "Running Simulation 64 MAC...");
        t64 = TileMatrix_GEMM_Num_Sim<64>(seqtile_A, seqtile_B, seqtile_c);
    }
    echo(stop_status, "");

    // printf("%d,%d,%d,%d,%d,%d,%lu,%lu,%lu,", seqtile->meta_m, seqtile->meta_n, seqtile->meta_nnz, seqtile_c->_nnz, seqtile_c->tile_val_map[seqtile_c->_nnz], seqtile_c->meta_nnz, t32.mp16, t32.mp4, t32.mpnz);
    printf("%d,%d,%d,%lu,%lu,%lu,", seqtile_A->meta_m, seqtile_A->meta_n, seqtile_A->meta_nnz, t32.mp16, t32.mp4, t32.mpnz);

    printf(
        "%lu,%.2lf,%lu,%.2lf,%lu,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%lu,%.2lf,%lu,%.2lf,%lu,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf\n",
        t32.ds_cycle, t32.ds_mac, t32.rm_cycle, t32.rm_mac, t32.smu_cycle, t32.smu_mac, 
        t32.ds_energy, t32.rm_energy, t32.smu_energy, t32.smu_speedup, 
        t64.ds_cycle, t64.ds_mac, t64.rm_cycle, t64.rm_mac, t64.smu_cycle, t64.smu_mac,
        t64.ds_energy, t64.rm_energy, t64.smu_energy, t64.smu_speedup,
        t32.ds_mt125, t32.ds_mt25, t32.ds_mt50, t32.ds_mt75, t32.ds_mt100, t32.rm_mt125, t32.rm_mt25, t32.rm_mt50, t32.rm_mt75, t32.rm_mt100, t32.smu_mt125, t32.smu_mt25, t32.smu_mt50, t32.smu_mt75, t32.smu_mt100,
        t64.ds_mt125, t64.ds_mt25, t64.ds_mt50, t64.ds_mt75, t64.ds_mt100, t64.rm_mt125, t64.rm_mt25, t64.rm_mt50, t64.rm_mt75, t64.rm_mt100, t64.smu_mt125, t64.smu_mt25, t64.smu_mt50, t64.smu_mt75, t64.smu_mt100,
        t32.ds_a_ntr, t32.ds_b_ntr, t32.ds_c_ntr, t32.rm_a_ntr, t32.rm_b_ntr, t32.rm_c_ntr, t32.smu_a_ntr, t32.smu_b_ntr, t32.smu_c_ntr, t64.ds_a_ntr, t64.ds_b_ntr, t64.ds_c_ntr, t64.rm_a_ntr, t64.rm_b_ntr, t64.rm_c_ntr, t64.smu_a_ntr, t64.smu_b_ntr, t64.smu_c_ntr,
        t32.ds_nt125, t32.ds_nt25, t32.ds_nt50, t32.ds_nt75, t32.ds_nt100, t32.rm_nt125, t32.rm_nt25, t32.rm_nt50, t32.rm_nt75, t32.rm_nt100, t32.smu_nt125, t32.smu_nt25, t32.smu_nt50, t32.smu_nt75, t32.smu_nt100,
        t64.ds_nt125, t64.ds_nt25, t64.ds_nt50, t64.ds_nt75, t64.ds_nt100, t64.rm_nt125, t64.rm_nt25, t64.rm_nt50, t64.rm_nt75, t64.rm_nt100, t64.smu_nt125, t64.smu_nt25, t64.smu_nt50, t64.smu_nt75, t64.smu_nt100,
        t32.rm_c_nut, t32.smu_a_nut, t32.smu_b_nut, t32.smu_c_nut, 
        t64.rm_c_nut, t64.smu_a_nut, t64.smu_b_nut, t64.smu_c_nut
    );
    fflush(stdout);
    destroy_TileMatrix(seqtile_A);
    if (seqtile_B != seqtile_A) destroy_TileMatrix(seqtile_B);
    destroy_TileMatrix(seqtile_c, true, false);

    return 0;
}