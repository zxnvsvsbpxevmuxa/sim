#pragma once
#include <tmatrix/Utils/timer.h>
#include <tmatrix/MMIO/mmio_highlevel.h>
#include <tmatrix/DataStructure/Factory.h>
#include <tmatrix/DataStructure/transform.h>
#include <tmatrix/Calculation/TileMatrix_x_v.h>

void test_SpMV(const char *filepath)
{
    MatIndex m, n, nnz;
    MatIndex *row_ptr, *col_ptr;
    MatValue *val_ptr;
    Timer read_mtx, csr_to_seqtile;
    bool exists;

    echo(title, "test SpMV on %s", filepath);

    char tile_path[128];
    sprintf(tile_path, "%s.tile", filepath);
    exists = file_exists(tile_path);
    if (!exists) {
        int isSymmetric;
        echo(start_status, "Reading Matrix file: [bold green]%s[/]", filepath);
        timer_start(read_mtx);
        mmio_allinone(filepath, &m, &n, &nnz, &isSymmetric, &row_ptr, &col_ptr, &val_ptr);
        timer_end(read_mtx);
        echo(stop_status, "");
        echo(success, "m = %d, n = %d, nnz = %d; time: %f ms", m, n, nnz, timer_duration(read_mtx));
    }
    echo(start_status, "CSR to TileMatrix");

    timer_start(csr_to_seqtile);
    TileMatrix *seqtile = CSR_to_TileMatrix(m, n, nnz, row_ptr, col_ptr, val_ptr, filepath);
    timer_end(csr_to_seqtile);
    echo(stop_status, "");
    echo(success, "CSR to TileMatrix time: %f ms", timer_duration(csr_to_seqtile));

    if (!exists) {
        free(row_ptr);
        free(col_ptr);
        free(val_ptr);
    }

    echo(start_status, "Prepare SpMV test data");
    MatValue *v = (MatValue *)malloc(n * sizeof(MatValue));
    MatValue *tile_res = (MatValue *)malloc(m * sizeof(MatValue));

#pragma omp parallel for
    for (MatIndex i = 0; i < n; ++i)
        v[i] = 1;
    echo(stop_status, "");
    echo(success, "SpMV test data prepared");
    
    echo(start_status, "Simulate SpMV 128 MAC");
    task_results t32 = TileMatrixSpMV<32>(seqtile);
    echo(start_status, "Simulate SpMV 64 MAC");
    task_results t64 = TileMatrixSpMV<64>(seqtile);
    echo(stop_status, "");
    
    free(v);
    free(tile_res);
    
    // print seqtile's meta info
    printf("%d,%d,%d,%lu,%lu,%lu,", seqtile->meta_m, seqtile->meta_n, seqtile->meta_nnz, t32.mp16, t32.mp4, t32.mpnz);

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
    destroy_TileMatrix(seqtile);
}
