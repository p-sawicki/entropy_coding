#ifndef ENTROPY_CODER_LIB_LOG_H
#define ENTROPY_CODER_LIB_LOG_H

#include <fstream>

namespace EntropyCoding {

enum class SyntaxElement {
  end_of_slice_one_bit = 0x00,
  end_of_tile_one_bit = 0x01,
  end_of_subset_one_bit = 0x02,

  alf_ctb_flag = 0x10,
  alf_use_aps_flag = 0x11,
  alf_luma_fixed_filter_idx = 0x12,
  alf_luma_prev_filter_idx = 0x13,
  alf_ctb_filter_alt_idx = 0x14,
  alf_ctb_cc_cb_idc = 0x15,
  alf_ctb_cc_cr_idx = 0x16,

  sao_merge_left_flag = 0x20,
  sao_merge_up_flag = 0x21,
  sao_type_idx_luma = 0x22,
  sao_type_idx_chroma = 0x23,
  sao_offset_abs = 0x24,
  sao_offset_sign_flag = 0x25,
  sao_band_position = 0x26,
  sao_eo_class_luma = 0x27,
  sao_eo_class_chroma = 0x28,

  split_cu_flag = 0x30,
  split_qt_flag = 0x31,
  mtt_split_cu_vertical_flag = 0x32,
  mtt_split_cu_binary_flag = 0x33,
  non_inter_flag = 0x34,

  cu_skip_flag = 0x40,
  pred_mode_ibc_flag = 0x41,
  pred_mode_plt_flag = 0x42,
  cu_act_enabled_flag = 0x43,
  pred_mode_flag = 0x44,
  intra_bdpcm_luma_flag = 0x45,
  intra_bdpcm_luma_dir_flag = 0x46,
  intra_mip_flag = 0x47,
  intra_mip_transposed_flag = 0x48,
  intra_mip_mode = 0x49,
  intra_luma_ref_idx = 0x4A,
  intra_subpartitions_mode_flag = 0x4B,
  intra_subpartitions_split_flag = 0x4C,
  intra_luma_mpm_flag = 0x4D,
  intra_luma_not_planar_flag = 0x4E,
  intra_luma_mpm_idx = 0x4F,
  intra_luma_mpm_remainder = 0x50,
  intra_bdpcm_chroma_flag = 0x51,
  intra_bdpcm_chroma_dir_flag = 0x52,
  cclm_mode_flag = 0x53,
  cclm_mode_idx = 0x54,
  intra_chroma_pred_mode = 0x55,
  general_merge_flag = 0x56,
  inter_pred_idc = 0x57,
  inter_affine_flag = 0x58,
  cu_affine_type_flag = 0x59,
  sym_mvd_flag = 0x5A,
  ref_idx_l0 = 0x5B,
  mvp_l0_flag = 0x5C,
  ref_idx_l1 = 0x5D,
  mvp_l1_flag = 0x5E,
  amvr_flag = 0x5F,
  amvr_precision_idx = 0x60,
  bcw_idx = 0x61,
  cu_coded_flag = 0x62,
  cu_sbt_flag = 0x63,
  cu_sbt_quad_flag = 0x64,
  cu_sbt_horizontal_flag = 0x65,
  cu_sbt_pos_flag = 0x66,
  lfnst_idx = 0x67,
  mts_idx = 0x68,

  palette_predictor_run = 0x70,
  num_signalled_palette_entries = 0x71,
  new_palette_entries = 0x72,
  palette_escape_val_present_flag = 0x73,
  palette_idx_idc = 0x74,
  palette_transpose_flag = 0x75,
  copy_above_palette_indices_flag = 0x76,
  run_copy_flag = 0x77,
  palette_escape_val = 0x78,

  regular_merge_flag = 0x80,
  mmvd_merge_flag = 0x81,
  mmvd_cand_flag = 0x82,
  mmvd_distance_idx = 0x83,
  mmvd_direction_idx = 0x84,
  ciip_flag = 0x85,
  merge_subblock_flag = 0x86,
  merge_subblock_idx = 0x87,
  merge_gpm_partition_idx = 0x88,
  merge_gpm_idx0 = 0x89,
  merge_gpm_idx1 = 0x8A,
  merge_idx = 0x8B,

  abs_mvd_greater0_flag = 0x90,
  abs_mvd_greater1_flag = 0x91,
  abs_mvd_minus2 = 0x92,
  mvd_sign_flag = 0x93,

  tu_y_coded_flag = 0xA0,
  tu_cb_coded_flag = 0xA1,
  tu_cr_coded_flag = 0xA2,
  cu_qp_delta_abs = 0xA3,
  cu_qp_delta_sign_flag = 0xA4,
  cu_chroma_qp_offset_flag = 0xA5,
  cu_chroma_qp_offset_idx = 0xA6,
  transform_skip_flag = 0xA7,
  tu_joint_cbcr_residual_flag = 0xA8,

  last_sig_coeff_x_prefix = 0xB0,
  last_sig_coeff_y_prefix = 0xB1,
  last_sig_coeff_x_suffix = 0xB2,
  last_sig_coeff_y_suffix = 0xB3,
  sb_coded_flag = 0xB4,
  sig_coeff_flag = 0xB5,
  par_level_flag = 0xB6,
  abs_level_gtx_flag = 0xB7,
  abs_remainder = 0xB8,
  dec_abs_level = 0xB9,
  coeff_sign_flag = 0xBA,
};

class Logger {
  ::std::ofstream fs;

  template <typename... vals, typename val>
  void LogValues(val first, vals... rest) {
    fs << static_cast<int>(first) << "\t";
    LogValues(rest...);
  }

  template <typename val> void LogValues(val first) {
    fs << static_cast<int>(first) << "\n";
  }

public:
  Logger(const char *filename) : fs(filename) {}

#ifdef ENABLE_LOGGING
  void LogElement(const SyntaxElement elem) {
    fs << static_cast<int>(elem) << "\n";
  }

  template <typename... vals>
  void LogElements(const SyntaxElement elem, vals... rest) {
    fs << static_cast<int>(elem) << "\t";
    LogValues(rest...);
  }

  void LogBits(const uint32_t bits) { fs << ::std::hex << bits; }

#else
  void LogElement(const SyntaxElement) {}
  template <typename... vals> void LogElements(const SyntaxElement, vals...) {}
  void LogBits(const uint32_t) {}
#endif // ENABLE_LOGGING

  Logger(const Logger &rhs) = delete;
  Logger &operator=(const Logger &rhs) = delete;
};

extern Logger binLogger, bitLogger;
} // namespace EntropyCoding

#endif // ENTROPY_CODER_LIB_LOG_H