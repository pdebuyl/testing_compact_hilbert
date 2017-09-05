/*
 * Test Hilbert space-filling curve for domain with unequal side lengths.
 * Copyright © 2013–2015 Peter Colberg.
 * Distributed under the MIT license. (See accompanying file LICENSE.)
 */

#include <Hilbert.hpp>
#include <hdf5.h>
#include <stdlib.h>

/**
Generates test data using [libhilbert] by Chris Hamilton.

Build libhilbert using `./configure` and `make`, then compile this program:

~~~
c++ -Ilibhilbert-0.2-1/include -O2 -o hilbert hilbert.cxx libhilbert-0.2-1/src/.libs/libHilbert.a -lhdf5
~~~

[libhilbert]: https://web.cs.dal.ca/~chamilto/hilbert/
**/
int main(int argc, char **argv)
{
  if (argc < 4 || argc > 6) {
    fprintf(stderr, "Usage: %s OUTPUT M1 M2 [M3 [M4]]\n", argv[0]);
    return 1;
  }
  int D = argc-2;
  int M[4];
  hsize_t L[4];
  hsize_t dims[2] = {1, D};
  for (int i = 0; i < 4; ++i) {
    M[i] = i < D ? atoi(argv[i+2]) : 0;
    L[i] = 1<<M[i];
    dims[0] *= L[i];
  }
  int *buf = (int *)malloc(dims[0]*dims[1]*sizeof(int));
  for (hsize_t x = 0; x < L[0]; ++x) {
    for (hsize_t y = 0; y < L[1]; ++y) {
      for (hsize_t z = 0; z < L[2]; ++z) {
        for (hsize_t w = 0; w < L[3]; ++w) {
          CFixBitVec p[4];
          p[0] = x;
          p[1] = y;
          p[2] = z;
          p[3] = w;
          CFixBitVec hc;
          Hilbert::coordsToCompactIndex(p, M, D, hc);
          int h = hc.rack();
          buf[h*D] = x;
          buf[h*D+1] = y;
          if (D >= 3) buf[h*D+2] = z;
          if (D >= 4) buf[h*D+3] = w;
        }
      }
    }
  }
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_libver_bounds(fapl, H5F_LIBVER_18, H5F_LIBVER_LATEST);
  hid_t file = H5Fcreate(argv[1], H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  hid_t space = H5Screate_simple(2, dims, NULL);
  hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl, 2, dims);
  H5Pset_shuffle(dcpl);
  H5Pset_deflate(dcpl, 9);
  hid_t dset = H5Dcreate(file, "value", H5T_NATIVE_UINT, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  H5Dwrite(dset, H5T_NATIVE_INT, space, H5S_ALL, H5P_DEFAULT, buf);
  hsize_t adims[1] = {D};
  H5Sset_extent_simple(space, 1, adims, NULL);
  hid_t attr = H5Acreate(dset, "M", H5T_NATIVE_UINT8, space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_INT, M);
  H5Pclose(dcpl);
  H5Pclose(fapl);
  H5Sclose(space);
  H5Dclose(dset);
  H5Aclose(attr);
  H5Fclose(file);
}
