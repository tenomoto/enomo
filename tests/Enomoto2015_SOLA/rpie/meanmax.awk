#1   87.722565515025522        1.8765646411143089E-014   7.7650950666798984E-016   4.0190073491430667E-014   2.2204460492503131E-016   8.0889057664263906E-033
{
  n++
  for (i = 3; i <=7; i++) { 
    mean[i] = mean[i] + $i
    max[i] = $i > max[i] ? $i : max[i]
  }
}
END {
  split(FILENAME, f, ".")
  ntrunc=substr(f[1], 8)
  printf("%d ", ntrunc)
  for (i = 3; i <=7; i++) { 
    printf("%e %e ", mean[i]/n, max[i])
  }
  printf("\n")
}
