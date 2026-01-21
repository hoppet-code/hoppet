#define DEFINE_SPLIT_FN(NAME) \
extern "C" { \
  double hoppet_cxx__splitfn__##NAME(double y, int piece); \
} \
namespace hoppet::split_fns { \
  inline double NAME(double y, int piece) { \
    return hoppet_cxx__splitfn__##NAME(y, piece); \
  } \
}

DEFINE_SPLIT_FN(pqq)
DEFINE_SPLIT_FN(pqg)
DEFINE_SPLIT_FN(pgq)
DEFINE_SPLIT_FN(pgg)

DEFINE_SPLIT_FN(p1qqv)
DEFINE_SPLIT_FN(p1qqbarv)
DEFINE_SPLIT_FN(p1qqs)
DEFINE_SPLIT_FN(p1qg)
DEFINE_SPLIT_FN(p1gq)
DEFINE_SPLIT_FN(p1gg)

DEFINE_SPLIT_FN(p2gg)
DEFINE_SPLIT_FN(p2qg2nf)
DEFINE_SPLIT_FN(p2ps)
DEFINE_SPLIT_FN(p2gq)
DEFINE_SPLIT_FN(p2nsplus)
DEFINE_SPLIT_FN(p2nsminus)
DEFINE_SPLIT_FN(p2nss)

DEFINE_SPLIT_FN(p3gg)
DEFINE_SPLIT_FN(p3qg2nf)
DEFINE_SPLIT_FN(p3ps)
DEFINE_SPLIT_FN(p3gq)
DEFINE_SPLIT_FN(p3nsplus)
DEFINE_SPLIT_FN(p3nsminus)
DEFINE_SPLIT_FN(p3nss)

#undef DEFINE_SPLIT_FN