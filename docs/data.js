/* Precomputed results. Every number here came out of the pipeline in
   homology_pipeline.py / the solubility notebooks; nothing is estimated. */

const BANDS = ["no hit", "<30%", "30-40%", "40-50%", "50-70%", "70-90%", ">90%"];

/* bands[dataset][split][features] = [{band, n, auc, acc, mcc}, ...]
   null  = that combination was not computed for this dataset. */
const RESULTS = {
  solubility: {
    label: "Solubility (DeepSol, 69,420 proteins)",
    random: {
      kmer: [
        { band: "no hit", n: 2330, auc: 0.611054, acc: 0.643777, mcc: 0.167130 },
        { band: "<30%",   n: 1259, auc: 0.684998, acc: 0.687847, mcc: 0.244597 },
        { band: "30-40%", n: 2766, auc: 0.687367, acc: 0.650759, mcc: 0.258618 },
        { band: "40-50%", n: 2034, auc: 0.702674, acc: 0.660275, mcc: 0.280810 },
        { band: "50-70%", n: 2479, auc: 0.687376, acc: 0.642194, mcc: 0.252963 },
        { band: "70-90%", n: 1526, auc: 0.730060, acc: 0.663172, mcc: 0.331341 },
        { band: ">90%",   n: 1491, auc: 0.517276, acc: 0.519115, mcc: 0.031843 }
      ],
      lookup: null,
      comp: [
        { band: "no hit", n: 2330, auc: 0.584507, acc: 0.633047, mcc: 0.134063 },
        { band: "<30%",   n: 1259, auc: 0.638652, acc: 0.675139, mcc: 0.184468 },
        { band: "30-40%", n: 2766, auc: 0.626449, acc: 0.604121, mcc: 0.137408 },
        { band: "40-50%", n: 2034, auc: 0.632383, acc: 0.620944, mcc: 0.183560 },
        { band: "50-70%", n: 2479, auc: 0.637197, acc: 0.603873, mcc: 0.159424 },
        { band: "70-90%", n: 1526, auc: 0.656690, acc: 0.612058, mcc: 0.231355 },
        { band: ">90%",   n: 1491, auc: 0.562713, acc: 0.531858, mcc: 0.064361 }
      ]
    },
    clustered: {
      kmer: [
        { band: "no hit", n: 3887, auc: 0.631665, acc: 0.641883, mcc: 0.208640 },
        { band: "<30%",   n: 3206, auc: 0.696412, acc: 0.678415, mcc: 0.269265 },
        { band: "30-40%", n: 3601, auc: 0.678874, acc: 0.633991, mcc: 0.241746 },
        { band: "40-50%", n: 1399, auc: 0.666971, acc: 0.626876, mcc: 0.231795 },
        { band: "50-70%", n:  961, auc: 0.687406, acc: 0.646202, mcc: 0.258624 },
        { band: "70-90%", n:  386, auc: 0.701804, acc: 0.639896, mcc: 0.296007 },
        { band: ">90%",   n:  445, auc: 0.521374, acc: 0.505618, mcc: 0.044224 }
      ],
      lookup: null,
      comp: [
        { band: "no hit", n: 3887, auc: 0.596392, acc: 0.613326, mcc: 0.132532 },
        { band: "<30%",   n: 3206, auc: 0.653081, acc: 0.646912, mcc: 0.169224 },
        { band: "30-40%", n: 3601, auc: 0.612827, acc: 0.590392, mcc: 0.136871 },
        { band: "40-50%", n: 1399, auc: 0.616376, acc: 0.592566, mcc: 0.152483 },
        { band: "50-70%", n:  961, auc: 0.613835, acc: 0.591051, mcc: 0.125254 },
        { band: "70-90%", n:  386, auc: 0.628084, acc: 0.544041, mcc: 0.107198 },
        { band: ">90%",   n:  445, auc: 0.513068, acc: 0.460674, mcc: -0.033773 }
      ]
    }
  },

  enzyme: {
    label: "Enzyme vs non-enzyme (UniProt, 101,304 proteins)",
    random: {
      kmer: [
        { band: "no hit", n: 1410, auc: 0.8143, acc: 0.7574, mcc: 0.3567 },
        { band: "<30%", n: 894, auc: 0.7688, acc: 0.7058, mcc: 0.4055 },
        { band: "30-40%", n: 2046, auc: 0.808, acc: 0.7493, mcc: 0.4989 },
        { band: "40-50%", n: 2023, auc: 0.8334, acc: 0.7533, mcc: 0.5076 },
        { band: "50-70%", n: 3851, auc: 0.8692, acc: 0.7936, mcc: 0.5885 },
        { band: "70-90%", n: 4899, auc: 0.9185, acc: 0.8392, mcc: 0.6764 },
        { band: ">90%", n: 5138, auc: 0.9436, acc: 0.8793, mcc: 0.7422 }
      ],
      comp: [
        { band: "no hit", n: 1410, auc: 0.7848, acc: 0.7355, mcc: 0.3223 },
        { band: "<30%", n: 894, auc: 0.7255, acc: 0.6846, mcc: 0.3572 },
        { band: "30-40%", n: 2046, auc: 0.7403, acc: 0.6711, mcc: 0.3412 },
        { band: "40-50%", n: 2023, auc: 0.7484, acc: 0.6733, mcc: 0.3527 },
        { band: "50-70%", n: 3851, auc: 0.7595, acc: 0.6726, mcc: 0.355 },
        { band: "70-90%", n: 4899, auc: 0.7718, acc: 0.6903, mcc: 0.3751 },
        { band: ">90%", n: 5138, auc: 0.759, acc: 0.7024, mcc: 0.3642 }
      ],
      lookup: [
        { band: "no hit", n: 1410, auc: null, acc: 0.8539, mcc: 0.0 },
        { band: "<30%", n: 894, auc: null, acc: 0.943, mcc: 0.8827 },
        { band: "30-40%", n: 2046, auc: null, acc: 0.957, mcc: 0.9147 },
        { band: "40-50%", n: 2023, auc: null, acc: 0.9733, mcc: 0.9464 },
        { band: "50-70%", n: 3851, auc: null, acc: 0.9914, mcc: 0.9827 },
        { band: "70-90%", n: 4899, auc: null, acc: 0.9957, mcc: 0.9914 },
        { band: ">90%", n: 5138, auc: null, acc: 0.9975, mcc: 0.9946 }
      ]
    },
    clustered: {
      kmer: [
        { band: "no hit", n: 4201, auc: 0.8348, acc: 0.7598, mcc: 0.4613 },
        { band: "<30%", n: 4258, auc: 0.8093, acc: 0.7428, mcc: 0.4876 },
        { band: "30-40%", n: 5574, auc: 0.8316, acc: 0.7569, mcc: 0.5157 },
        { band: "40-50%", n: 2839, auc: 0.826, acc: 0.7372, mcc: 0.4785 },
        { band: "50-70%", n: 2445, auc: 0.7651, acc: 0.7014, mcc: 0.3815 },
        { band: "70-90%", n: 718, auc: 0.88, acc: 0.7813, mcc: 0.5504 },
        { band: ">90%", n: 226, auc: 0.8228, acc: 0.7345, mcc: 0.3672 }
      ],
      comp: [
        { band: "no hit", n: 4201, auc: 0.786, acc: 0.7222, mcc: 0.3842 },
        { band: "<30%", n: 4258, auc: 0.7181, acc: 0.6768, mcc: 0.3558 },
        { band: "30-40%", n: 5574, auc: 0.7539, acc: 0.676, mcc: 0.3608 },
        { band: "40-50%", n: 2839, auc: 0.7658, acc: 0.6675, mcc: 0.3414 },
        { band: "50-70%", n: 2445, auc: 0.6905, acc: 0.6335, mcc: 0.2345 },
        { band: "70-90%", n: 718, auc: 0.797, acc: 0.688, mcc: 0.3513 },
        { band: ">90%", n: 226, auc: 0.7384, acc: 0.6504, mcc: 0.1649 }
      ],
      lookup: [
        { band: "no hit", n: 4201, auc: null, acc: 0.7532, mcc: 0.0 },
        { band: "<30%", n: 4258, auc: null, acc: 0.8953, mcc: 0.7979 },
        { band: "30-40%", n: 5574, auc: null, acc: 0.9376, mcc: 0.8744 },
        { band: "40-50%", n: 2839, auc: null, acc: 0.9616, mcc: 0.9238 },
        { band: "50-70%", n: 2445, auc: null, acc: 0.9779, mcc: 0.9535 },
        { band: "70-90%", n: 718, auc: null, acc: 0.9916, mcc: 0.9829 },
        { band: ">90%", n: 226, auc: null, acc: 0.9646, mcc: 0.9215 }
      ]
    }
  }
};

/* Fraction of best-hit pairs where the test protein and its nearest
   training protein carry the SAME label. The key diagnostic. */
const AGREEMENT = {
  solubility: {
    random: [
      { band: "<30%",   n: 1259, v: 0.617156 },
      { band: "30-40%", n: 2766, v: 0.652205 },
      { band: "40-50%", n: 2034, v: 0.681416 },
      { band: "50-70%", n: 2479, v: 0.713998 },
      { band: "70-90%", n: 1526, v: 0.749017 },
      { band: ">90%",   n: 1491, v: 0.099262 }
    ],
    clustered: null
  },
  enzyme: {
    random: [
      { band: "<30%",   n:  894, v: 0.943 },
      { band: "30-40%", n: 2046, v: 0.957 },
      { band: "40-50%", n: 2023, v: 0.973 },
      { band: "50-70%", n: 3851, v: 0.991 },
      { band: "70-90%", n: 4899, v: 0.996 },
      { band: ">90%",   n: 5138, v: 0.998 }
    ],
    clustered: [
      { band: "<30%",   n: 4258, v: 0.895 },
      { band: "30-40%", n: 5574, v: 0.938 },
      { band: "40-50%", n: 2839, v: 0.962 },
      { band: "50-70%", n: 2445, v: 0.978 },
      { band: "70-90%", n:  718, v: 0.992 },
      { band: ">90%",   n:  226, v: 0.965 }
    ]
  }
};
