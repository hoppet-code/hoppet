#define DEFINE_SPLIT_FN(NAME) \
extern "C" { \
  double hoppet_cxx__split_fn__##NAME(double y, int piece); \
} \
namespace hoppet::split_fns { \
  inline double NAME(double y, int piece) { \
    return hoppet_cxx__split_fn__##NAME(y, piece); \
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


#undef DEFINE_SPLIT_FN