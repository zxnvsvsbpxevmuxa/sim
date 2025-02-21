#include <msg.h>

#ifdef SPMV
#include <test/SpMV.h>
int main(int argc, char **argv)
{
    test_SpMV(argv[1]);
    return 0;
}
#endif

#ifdef SPGEMM
#include <test/SpGEMM.h>
int main(int argc, char **argv)
{
    return test_SpGEMM(argv[1]);
}
#endif

#ifdef SPGEMM2
#include <test/SpGEMM.h>
int main(int argc, char **argv)
{
    return test_SpGEMM(argv[1], argv[2]);
}
#endif

#ifdef SPMSPV
#include <test/SpMSpV.h>
int main(int argc, char **argv)
{
    test_SpMSpV(argv[1]);
    return 0;
}
#endif

#ifdef SPMM
#include <test/SpMM.h>
int main(int argc, char **argv)
{
    test_SpMM(argv[1]);
    return 0;
}
#endif
