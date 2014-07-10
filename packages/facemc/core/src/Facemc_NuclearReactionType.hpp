//---------------------------------------------------------------------------//
//! 
//! \file   Facemc_NuclearReactionType.hpp
//! \author Alex Robinson
//! \brief  Nuclear reaction type enumeration and helper function declarations.
//!
//---------------------------------------------------------------------------//

#ifndef FACEMC_NUCLEAR_REACTION_TYPE_HPP
#define FACEMC_NUCLEAR_REACTION_TYPE_HPP

namespace Facemc{

/*! The nuclear reaction type enum.
 * \details Converting the enumeration value to an integer will recover the
 * corresponding MT #.
 */
enum NuclearReactionType{
  N__TOTAL_REACTION = 1,
  N__N_ELASTIC_REACTION = 2,
  N__N_NON_ELASTIC_REACTION = 3, // sum of 4,16,17,22-37,41,42,44,45,152-154,156-181,183-190,194-196,198-200
  N__N_INELASTIC_REACTION = 4, // sum of 50-91
  N__ANYTHING_REACTION = 5,
  N__2N_D_REACTION = 11,
  N__2N_REACTION = 16,
  N__3N_REACTION = 17,
  N__TOTAL_FISSION_REACTION = 18, // sum of 19-21,38
  N__FISSION_REACTION = 19,
  N__N_FISSION_REACTION = 20,
  N__2N_FISSION_REACTION = 21,
  N__N_ALPHA_REACTION = 22,
  N__N_3ALPHA_REACTION = 23,
  N__2N_ALPHA_REACTION = 24,
  N__3N_ALPHA_REACTION = 25,
  N__TOTAL_ABSORPTION = 27, // sum 0f 18, 101
  N__N_P_REACTION = 28,
  N__N_2ALPHA_REACTION = 29,
  N__2N_2ALPHA_REACTION = 30,
  N__N_D_REACTION = 32,
  N__N_T_REACTION = 33,
  N__N_HE3_REACTION = 34,
  N__N_D_2ALPHA_REACTION = 35,
  N__N_T_2ALPHA_REACTION = 36,
  N__4N_REACTION = 37,
  N__3N_FISSION_REACTION = 38,
  N__2N_P_REACTION = 41,
  N__3N_P_REACTION = 42,
  N__N_2P_REACTION = 44,
  N__N_P_ALPHA_REACTION = 45,
  N__N_EXICTED_STATE_1_REACTION = 51,
  N__N_EXCITED_STATE_2_REACTION = 52,
  N__N_EXCITED_STATE_3_REACTION = 53,
  N__N_EXCITED_STATE_4_REACTION = 54,
  N__N_EXCITED_STATE_5_REACTION = 55,
  N__N_EXCITED_STATE_6_REACTION = 56,
  N__N_EXCITED_STATE_7_REACTION = 57,
  N__N_EXCITED_STATE_8_REACTION = 58,
  N__N_EXCITED_STATE_9_REACTION = 59,
  N__N_EXCITED_STATE_10_REACTION = 60,
  N__N_EXCITED_STATE_11_REACTION = 61,
  N__N_EXCITED_STATE_12_REACTION = 62,
  N__N_EXCITED_STATE_13_REACTION = 63,
  N__N_EXCITED_STATE_14_REACTION = 64,
  N__N_EXCITED_STATE_15_REACTION = 65,
  N__N_EXCITED_STATE_16_REACTION = 66,
  N__N_EXCITED_STATE_17_REACTION = 67,
  N__N_EXCITED_STATE_18_REACTION = 68,
  N__N_EXCITED_STATE_19_REACTION = 69,
  N__N_EXCITED_STATE_20_REACTION = 70,
  N__N_EXCITED_STATE_21_REACTION = 71,
  N__N_EXCITED_STATE_22_REACTION = 72,
  N__N_EXCITED_STATE_23_REACTION = 73,
  N__N_EXCITED_STATE_24_REACTION = 74,
  N__N_EXCITED_STATE_25_REACTION = 75,
  N__N_EXCITED_STATE_26_REACTION = 76,
  N__N_EXCITED_STATE_27_REACTION = 77,
  N__N_EXCITED_STATE_28_REACTION = 78,
  N__N_EXCITED_STATE_29_REACTION = 79,
  N__N_EXCITED_STATE_30_REACTION = 80,
  N__N_EXCITED_STATE_31_REACTION = 81,
  N__N_EXCITED_STATE_32_REACTION = 82,
  N__N_EXCITED_STATE_33_REACTION = 83,
  N__N_EXCITED_STATE_34_REACTION = 84,
  N__N_EXCITED_STATE_35_REACTION = 85,
  N__N_EXCITED_STATE_36_REACTION = 86,
  N__N_EXCITED_STATE_37_REACTION = 87,
  N__N_EXCITED_STATE_38_REACTION = 88,
  N__N_EXCITED_STATE_39_REACTION = 89,
  N__N_EXCITED_STATE_40_REACTION = 90,
  N__N_CONTINUUM_REACTION = 91,
  N__CAPTURE_REACTION = 101, // sum of 102-117
  N__GAMMA_REACTION = 102,
  N__P_REACTION = 103,
  N__D_REACTION = 104,
  N__T_REACTION = 105,
  N__HE3_REACTION = 106,
  N__ALPHA_REACTION = 107,
  N__2ALPHA_REACTION = 108,
  N__3ALPHA_REACTION = 109,
  N__2P_REACTION = 111,
  N__P_ALPHA_REACTION = 112,
  N__T_2ALPHA_REACTION = 113,
  N__D_2ALPHA_REACTION = 114,
  N__P_D_REACTION = 115,
  N__P_T_REACTION = 116,
  N__D_ALPHA_REACTION = 117,
  N__5N_REACTION = 152,
  N__6N_REACTION = 153,
  N__2N_T_REACTION = 154,
  N__T_ALPHA_REACTION = 155,
  N__4N_P_REACTION = 156,
  N__3N_D_REACTION = 157,
  N__N_D_ALPHA_REACTION = 158,
  N__2N_P_ALPHA_REACTION = 159,
  N__7N_REACTION = 160,
  N__8N_REACTION = 161,
  N__5N_P_REACTION = 162,
  N__6N_P_REACTION = 163,
  N__7N_P_REACTION = 164,
  N__4N_ALPHA_REACTION = 165,
  N__5N_ALPHA_REACTION = 166,
  N__6N_ALPHA_REACTION = 167,
  N__7N_ALPHA_REACTION = 168,
  N__4N_D_REACTION = 169,
  N__5N_D_REACTION = 170,
  N__6N_D_REACTION = 171,
  N__3N_T_REACTION = 172,
  N__4N_T_REACTION = 173,
  N__5N_T_REACTION = 174,
  N__6N_T_REACTION = 175,
  N__2N_HE3_REACTION = 176,
  N__3N_HE3_REACTION = 177,
  N__4N_HE3_REACTION = 178,
  N__3N_2P_REACTION = 179,
  N__3N_2ALPHA_REACTION = 180,
  N__3N_P_ALPHA_REACTION = 181,
  N__D_T_REACTION = 182,
  N__N_P_D_REACTION = 183,
  N__N_P_T_REACTION = 184,
  N__N_D_T_REACTION = 185,
  N__N_P_HE3_REACTION = 186,
  N__N_D_HE3_REACTION = 187,
  N__N_T_HE3_REACTION = 188,
  N__N_T_ALPHA_REACTION = 189,
  N__2N_2P_REACTION = 190,
  N__P_HE3_REACTION = 191,
  N__D_HE3_REACTION = 192,
  N__HE3_ALPHA_REACTION = 193,
  N__4N_2P_REACTION = 194,
  N__4N_2ALPHA_REACTION = 195,
  N__4N_P_ALPHA_REACTION = 196,
  N__3P_REACTION = 197,
  N__N_3P_REACTION = 198,
  N__3N_2P_ALPHA_REACTION = 199,
  N__5N_2P_REACTION = 200,
  N__TOTAL_N_PRODUCTION = 201,
  N__TOTAL_GAMMA_PRODUCTION = 202,
  N__TOTAL_P_PRODUCTION = 203,
  N__TOTAL_D_PRODUCTION = 204,
  N__TOTAL_T_PRODUCTION = 205,
  N__TOTAL_HE3_PRODUCTION = 206,
  N__TOTAL_ALPHA_PRODUCTION = 207,
  N__AVERAGE_HEATING = 301,
  N__DPA = 444,
  N__P_EXCITED_STATE_0_REACTION = 600,
  N__P_EXCITED_STATE_1_REACTION = 601,
  N__P_EXCITED_STATE_2_REACTION = 602,
  N__P_EXCITED_STATE_3_REACTION = 603,
  N__P_EXCITED_STATE_4_REACTION = 604,
  N__P_EXCITED_STATE_5_REACTION = 605,
  N__P_EXCITED_STATE_6_REACTION = 606,
  N__P_EXCITED_STATE_7_REACTION = 607,
  N__P_EXCITED_STATE_8_REACTION = 608,
  N__P_EXCITED_STATE_9_REACTION = 609,
  N__P_EXCITED_STATE_10_REACTION = 610,
  N__P_EXCITED_STATE_11_REACTION = 611,
  N__P_EXCITED_STATE_12_REACTION = 612,
  N__P_EXCITED_STATE_13_REACTION = 613,
  N__P_EXCITED_STATE_14_REACTION = 614,
  N__P_EXCITED_STATE_15_REACTION = 615,
  N__P_EXCITED_STATE_16_REACTION = 616,
  N__P_EXCITED_STATE_17_REACTION = 617,
  N__P_EXCITED_STATE_18_REACTION = 618,
  N__P_EXCITED_STATE_19_REACTION = 619,
  N__P_EXCITED_STATE_20_REACTION = 620,
  N__P_EXCITED_STATE_21_REACTION = 621,
  N__P_EXCITED_STATE_22_REACTION = 622,
  N__P_EXCITED_STATE_23_REACTION = 623,
  N__P_EXCITED_STATE_24_REACTION = 624,
  N__P_EXCITED_STATE_25_REACTION = 625,
  N__P_EXCITED_STATE_26_REACTION = 626,
  N__P_EXCITED_STATE_27_REACTION = 627,
  N__P_EXCITED_STATE_28_REACTION = 628,
  N__P_EXCITED_STATE_29_REACTION = 629,
  N__P_EXCITED_STATE_30_REACTION = 630,
  N__P_EXCITED_STATE_31_REACTION = 631,
  N__P_EXCITED_STATE_32_REACTION = 632,
  N__P_EXCITED_STATE_33_REACTION = 633,
  N__P_EXCITED_STATE_34_REACTION = 634,
  N__P_EXCITED_STATE_35_REACTION = 635,
  N__P_EXCITED_STATE_36_REACTION = 636,
  N__P_EXCITED_STATE_37_REACTION = 637,
  N__P_EXCITED_STATE_38_REACTION = 638,
  N__P_EXCITED_STATE_39_REACTION = 639,
  N__P_EXCITED_STATE_40_REACTION = 640,
  N__P_EXCITED_STATE_41_REACTION = 641,
  N__P_EXCITED_STATE_42_REACTION = 642,
  N__P_EXCITED_STATE_43_REACTION = 643,
  N__P_EXCITED_STATE_44_REACTION = 644,
  N__P_EXCITED_STATE_45_REACTION = 645,
  N__P_EXCITED_STATE_46_REACTION = 646,
  N__P_EXCITED_STATE_47_REACTION = 647,
  N__P_EXCITED_STATE_48_REACTION = 648,
  N__P_CONTINUUM_REACTION = 649,
  N__D_EXCITED_STATE_0_REACTION = 650,
  N__D_EXCITED_STATE_1_REACTION = 651,
  N__D_EXCITED_STATE_2_REACTION = 652,
  N__D_EXCITED_STATE_3_REACTION = 653,
  N__D_EXCITED_STATE_4_REACTION = 654,
  N__D_EXCITED_STATE_5_REACTION = 655,
  N__D_EXCITED_STATE_6_REACTION = 656,
  N__D_EXCITED_STATE_7_REACTION = 657,
  N__D_EXCITED_STATE_8_REACTION = 658,
  N__D_EXCITED_STATE_9_REACTION = 659,
  N__D_EXCITED_STATE_10_REACTION = 660,
  N__D_EXCITED_STATE_11_REACTION = 661,
  N__D_EXCITED_STATE_12_REACTION = 662,
  N__D_EXCITED_STATE_13_REACTION = 663,
  N__D_EXCITED_STATE_14_REACTION = 664,
  N__D_EXCITED_STATE_15_REACTION = 665,
  N__D_EXCITED_STATE_16_REACTION = 666,
  N__D_EXCITED_STATE_17_REACTION = 667,
  N__D_EXCITED_STATE_18_REACTION = 668,
  N__D_EXCITED_STATE_19_REACTION = 669,
  N__D_EXCITED_STATE_20_REACTION = 670,
  N__D_EXCITED_STATE_21_REACTION = 671,
  N__D_EXCITED_STATE_22_REACTION = 672,
  N__D_EXCITED_STATE_23_REACTION = 673,
  N__D_EXCITED_STATE_24_REACTION = 674,
  N__D_EXCITED_STATE_25_REACTION = 675,
  N__D_EXCITED_STATE_26_REACTION = 676,
  N__D_EXCITED_STATE_27_REACTION = 677,
  N__D_EXCITED_STATE_28_REACTION = 678,
  N__D_EXCITED_STATE_29_REACTION = 679,
  N__D_EXCITED_STATE_30_REACTION = 680,
  N__D_EXCITED_STATE_31_REACTION = 681,
  N__D_EXCITED_STATE_32_REACTION = 682,
  N__D_EXCITED_STATE_33_REACTION = 683,
  N__D_EXCITED_STATE_34_REACTION = 684,
  N__D_EXCITED_STATE_35_REACTION = 685,
  N__D_EXCITED_STATE_36_REACTION = 686,
  N__D_EXCITED_STATE_37_REACTION = 687,
  N__D_EXCITED_STATE_38_REACTION = 688,
  N__D_EXCITED_STATE_39_REACTION = 689,
  N__D_EXCITED_STATE_40_REACTION = 690,
  N__D_EXCITED_STATE_41_REACTION = 691,
  N__D_EXCITED_STATE_42_REACTION = 692,
  N__D_EXCITED_STATE_43_REACTION = 693,
  N__D_EXCITED_STATE_44_REACTION = 694,
  N__D_EXCITED_STATE_45_REACTION = 695,
  N__D_EXCITED_STATE_46_REACTION = 696,
  N__D_EXCITED_STATE_47_REACTION = 697,
  N__D_EXCITED_STATE_48_REACTION = 698,
  N__D_CONTINUUM_REACTION = 699,
  N__T_EXCITED_STATE_0_REACTION = 700,
  N__T_EXCITED_STATE_1_REACTION = 701,
  N__T_EXCITED_STATE_2_REACTION = 702,
  N__T_EXCITED_STATE_3_REACTION = 703,
  N__T_EXCITED_STATE_4_REACTION = 704,
  N__T_EXCITED_STATE_5_REACTION = 705,
  N__T_EXCITED_STATE_6_REACTION = 706,
  N__T_EXCITED_STATE_7_REACTION = 707,
  N__T_EXCITED_STATE_8_REACTION = 708,
  N__T_EXCITED_STATE_9_REACTION = 709,
  N__T_EXCITED_STATE_10_REACTION = 710,
  N__T_EXCITED_STATE_11_REACTION = 711,
  N__T_EXCITED_STATE_12_REACTION = 712,
  N__T_EXCITED_STATE_13_REACTION = 713,
  N__T_EXCITED_STATE_14_REACTION = 714,
  N__T_EXCITED_STATE_15_REACTION = 715,
  N__T_EXCITED_STATE_16_REACTION = 716,
  N__T_EXCITED_STATE_17_REACTION = 717,
  N__T_EXCITED_STATE_18_REACTION = 718,
  N__T_EXCITED_STATE_19_REACTION = 719,
  N__T_EXCITED_STATE_20_REACTION = 720,
  N__T_EXCITED_STATE_21_REACTION = 721,
  N__T_EXCITED_STATE_22_REACTION = 722,
  N__T_EXCITED_STATE_23_REACTION = 723,
  N__T_EXCITED_STATE_24_REACTION = 724,
  N__T_EXCITED_STATE_25_REACTION = 725,
  N__T_EXCITED_STATE_26_REACTION = 726,
  N__T_EXCITED_STATE_27_REACTION = 727,
  N__T_EXCITED_STATE_28_REACTION = 728,
  N__T_EXCITED_STATE_29_REACTION = 729,
  N__T_EXCITED_STATE_30_REACTION = 730,
  N__T_EXCITED_STATE_31_REACTION = 731,
  N__T_EXCITED_STATE_32_REACTION = 732,
  N__T_EXCITED_STATE_33_REACTION = 733,
  N__T_EXCITED_STATE_34_REACTION = 734,
  N__T_EXCITED_STATE_35_REACTION = 735,
  N__T_EXCITED_STATE_36_REACTION = 736,
  N__T_EXCITED_STATE_37_REACTION = 737,
  N__T_EXCITED_STATE_38_REACTION = 738,
  N__T_EXCITED_STATE_39_REACTION = 739,
  N__T_EXCITED_STATE_40_REACTION = 740,
  N__T_EXCITED_STATE_41_REACTION = 741,
  N__T_EXCITED_STATE_42_REACTION = 742,
  N__T_EXCITED_STATE_43_REACTION = 743,
  N__T_EXCITED_STATE_44_REACTION = 744,
  N__T_EXCITED_STATE_45_REACTION = 745,
  N__T_EXCITED_STATE_46_REACTION = 746,
  N__T_EXCITED_STATE_47_REACTION = 747,
  N__T_EXCITED_STATE_48_REACTION = 748,
  N__T_CONTINUUM_REACTION = 749,
  N__HE3_EXCITED_STATE_0_REACTION = 750,
  N__HE3_EXCITED_STATE_1_REACTION = 751,
  N__HE3_EXCITED_STATE_2_REACTION = 752,
  N__HE3_EXCITED_STATE_3_REACTION = 753,
  N__HE3_EXCITED_STATE_4_REACTION = 754,
  N__HE3_EXCITED_STATE_5_REACTION = 755,
  N__HE3_EXCITED_STATE_6_REACTION = 756,
  N__HE3_EXCITED_STATE_7_REACTION = 757,
  N__HE3_EXCITED_STATE_8_REACTION = 758,
  N__HE3_EXCITED_STATE_9_REACTION = 759,
  N__HE3_EXCITED_STATE_10_REACTION = 760,
  N__HE3_EXCITED_STATE_11_REACTION = 761,
  N__HE3_EXCITED_STATE_12_REACTION = 762,
  N__HE3_EXCITED_STATE_13_REACTION = 763,
  N__HE3_EXCITED_STATE_14_REACTION = 764,
  N__HE3_EXCITED_STATE_15_REACTION = 765,
  N__HE3_EXCITED_STATE_16_REACTION = 766,
  N__HE3_EXCITED_STATE_17_REACTION = 767,
  N__HE3_EXCITED_STATE_18_REACTION = 768,
  N__HE3_EXCITED_STATE_19_REACTION = 769,
  N__HE3_EXCITED_STATE_20_REACTION = 770,
  N__HE3_EXCITED_STATE_21_REACTION = 771,
  N__HE3_EXCITED_STATE_22_REACTION = 772,
  N__HE3_EXCITED_STATE_23_REACTION = 773,
  N__HE3_EXCITED_STATE_24_REACTION = 774,
  N__HE3_EXCITED_STATE_25_REACTION = 775,
  N__HE3_EXCITED_STATE_26_REACTION = 776,
  N__HE3_EXCITED_STATE_27_REACTION = 777,
  N__HE3_EXCITED_STATE_28_REACTION = 778,
  N__HE3_EXCITED_STATE_29_REACTION = 779,
  N__HE3_EXCITED_STATE_30_REACTION = 780,
  N__HE3_EXCITED_STATE_31_REACTION = 781,
  N__HE3_EXCITED_STATE_32_REACTION = 782,
  N__HE3_EXCITED_STATE_33_REACTION = 783,
  N__HE3_EXCITED_STATE_34_REACTION = 784,
  N__HE3_EXCITED_STATE_35_REACTION = 785,
  N__HE3_EXCITED_STATE_36_REACTION = 786,
  N__HE3_EXCITED_STATE_37_REACTION = 787,
  N__HE3_EXCITED_STATE_38_REACTION = 788,
  N__HE3_EXCITED_STATE_39_REACTION = 789,
  N__HE3_EXCITED_STATE_40_REACTION = 790,
  N__HE3_EXCITED_STATE_41_REACTION = 791,
  N__HE3_EXCITED_STATE_42_REACTION = 792,
  N__HE3_EXCITED_STATE_43_REACTION = 793,
  N__HE3_EXCITED_STATE_44_REACTION = 794,
  N__HE3_EXCITED_STATE_45_REACTION = 795,
  N__HE3_EXCITED_STATE_46_REACTION = 796,
  N__HE3_EXCITED_STATE_47_REACTION = 797,
  N__HE3_EXCITED_STATE_48_REACTION = 798,
  N__HE3_CONTINUUM_REACTION = 799,
  N__ALPHA_EXCITED_STATE_0_REACTION = 800,
  N__ALPHA_EXCITED_STATE_1_REACTION = 801,
  N__ALPHA_EXCITED_STATE_2_REACTION = 802,
  N__ALPHA_EXCITED_STATE_3_REACTION = 803,
  N__ALPHA_EXCITED_STATE_4_REACTION = 804,
  N__ALPHA_EXCITED_STATE_5_REACTION = 805,
  N__ALPHA_EXCITED_STATE_6_REACTION = 806,
  N__ALPHA_EXCITED_STATE_7_REACTION = 807,
  N__ALPHA_EXCITED_STATE_8_REACTION = 808,
  N__ALPHA_EXCITED_STATE_9_REACTION = 809,
  N__ALPHA_EXCITED_STATE_10_REACTION = 810,
  N__ALPHA_EXCITED_STATE_11_REACTION = 811,
  N__ALPHA_EXCITED_STATE_12_REACTION = 812,
  N__ALPHA_EXCITED_STATE_13_REACTION = 813,
  N__ALPHA_EXCITED_STATE_14_REACTION = 814,
  N__ALPHA_EXCITED_STATE_15_REACTION = 815,
  N__ALPHA_EXCITED_STATE_16_REACTION = 816,
  N__ALPHA_EXCITED_STATE_17_REACTION = 817,
  N__ALPHA_EXCITED_STATE_18_REACTION = 818,
  N__ALPHA_EXCITED_STATE_19_REACTION = 819,
  N__ALPHA_EXCITED_STATE_20_REACTION = 820,
  N__ALPHA_EXCITED_STATE_21_REACTION = 821,
  N__ALPHA_EXCITED_STATE_22_REACTION = 822,
  N__ALPHA_EXCITED_STATE_23_REACTION = 823,
  N__ALPHA_EXCITED_STATE_24_REACTION = 824,
  N__ALPHA_EXCITED_STATE_25_REACTION = 825,
  N__ALPHA_EXCITED_STATE_26_REACTION = 826,
  N__ALPHA_EXCITED_STATE_27_REACTION = 827,
  N__ALPHA_EXCITED_STATE_28_REACTION = 828,
  N__ALPHA_EXCITED_STATE_29_REACTION = 829,
  N__ALPHA_EXCITED_STATE_30_REACTION = 830,
  N__ALPHA_EXCITED_STATE_31_REACTION = 831,
  N__ALPHA_EXCITED_STATE_32_REACTION = 832,
  N__ALPHA_EXCITED_STATE_33_REACTION = 833,
  N__ALPHA_EXCITED_STATE_34_REACTION = 834,
  N__ALPHA_EXCITED_STATE_35_REACTION = 835,
  N__ALPHA_EXCITED_STATE_36_REACTION = 836,
  N__ALPHA_EXCITED_STATE_37_REACTION = 837,
  N__ALPHA_EXCITED_STATE_38_REACTION = 838,
  N__ALPHA_EXCITED_STATE_39_REACTION = 839,
  N__ALPHA_EXCITED_STATE_40_REACTION = 840,
  N__ALPHA_EXCITED_STATE_41_REACTION = 841,
  N__ALPHA_EXCITED_STATE_42_REACTION = 842,
  N__ALPHA_EXCITED_STATE_43_REACTION = 843,
  N__ALPHA_EXCITED_STATE_44_REACTION = 844,
  N__ALPHA_EXCITED_STATE_45_REACTION = 845,
  N__ALPHA_EXCITED_STATE_46_REACTION = 846,
  N__ALPHA_EXCITED_STATE_47_REACTION = 847,
  N__ALPHA_EXCITED_STATE_48_REACTION = 848,
  N__ALPHA_CONTINUUM_REACTION = 849,
  N__2N_EXCITED_STATE_0_REACTION = 875,
  N__2N_EXCITED_STATE_1_REACTION = 876,
  N__2N_EXCITED_STATE_2_REACTION = 877,
  N__2N_EXCITED_STATE_3_REACTION = 878,
  N__2N_EXCITED_STATE_4_REACTION = 879,
  N__2N_EXCITED_STATE_5_REACTION = 880,
  N__2N_EXCITED_STATE_6_REACTION = 881,
  N__2N_EXCITED_STATE_7_REACTION = 882,
  N__2N_EXCITED_STATE_8_REACTION = 883,
  N__2N_EXCITED_STATE_9_REACTION = 884,
  N__2N_EXCITED_STATE_10_REACTION = 885,
  N__2N_EXCITED_STATE_11_REACTION = 886,
  N__2N_EXCITED_STATE_12_REACTION = 887,
  N__2N_EXCITED_STATE_13_REACTION = 888,
  N__2N_EXCITED_STATE_14_REACTION = 889,
  N__2N_EXCITED_STATE_15_REACTION = 890,
  N__2N_CONTINUUM_REACTION = 891,
};

//! Convert an unsigned int to a NuclearReactionType
NuclearReactionType convertUnsignedToNuclearReactionType( 
						     const unsigned reaction );

} // end Facemc namespace

#endif // end FACEMC_NUCLEAR_REACTION_TYPE_HPP

//---------------------------------------------------------------------------//
// end Facemc_NuclearReactionType.hpp
//---------------------------------------------------------------------------//
