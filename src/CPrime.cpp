#include <CPrime.h>

#include <cassert>
#include <cmath>
#include <set>

#define CPrimeMgrInst CPrimeMgr::instance()

class CPrimeMgr {
 public:
  static CPrimeMgr *instance();

  bool isPrime(int i) const;

  std::vector<int> factors(int n) const;

 private:
  CPrimeMgr() { }
 ~CPrimeMgr() { }

  void initPrimes();

  void genPrime(int n) const;

  bool checkPrime(int n) const;

 private:
  typedef std::set<int> Primes;

  mutable Primes primes_;
  mutable int    prime_max_;
};

CPrimeMgr *
CPrimeMgr::
instance()
{
  static CPrimeMgr *inst;

  if (! inst) {
    inst = new CPrimeMgr;

    inst->initPrimes();
  }

  return inst;
}

void
CPrimeMgr::
initPrimes()
{
  primes_.insert(2);
  primes_.insert(3);
  primes_.insert(5);
  primes_.insert(7);

  prime_max_ = 10;
}

bool
CPrimeMgr::
isPrime(int i) const
{
  assert(i > 0);

  if (i == 1) return true;

  genPrime(i);

  return (primes_.find(i) != primes_.end());
}

void
CPrimeMgr::
genPrime(int n) const
{
  if (n <= prime_max_)
    return;

  for (int i = prime_max_ + 1; i <= n; ++i) {
    if (checkPrime(i))
      primes_.insert(i);
  }

  prime_max_ = n;
}

bool
CPrimeMgr::
checkPrime(int n) const
{
  for (Primes::const_iterator p = primes_.begin(); p != primes_.end(); ++p) {
    int i = *p;

    if (n % i == 0) return false;
  }

  return true;
}

std::vector<int>
CPrimeMgr::
factors(int n) const
{
  std::vector<int> v;

  if (isPrime(n))
    v.push_back(n);
  else {
    int n1 = n/2;

    for (int i = n1; i >= 2; --i) {
      if (! isPrime(i)) continue;

      if (n % i != 0) continue;

      v.push_back(i);

      std::vector<int> v1 = factors(n / i);

      std::copy(v1.begin(), v1.end(), std::back_inserter(v));

      break;
    }
  }

  return v;
}

//---

bool
CPrime::
isPrime(int i)
{
  return CPrimeMgrInst->isPrime(i);
}

std::vector<int>
CPrime::
factors(int n)
{
  return CPrimeMgrInst->factors(n);
}
