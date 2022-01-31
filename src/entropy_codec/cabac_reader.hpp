#ifndef ENTROPY_CODEC_CABAC_READER
#define ENTROPY_CODEC_CABAC_READER

#include "arith_codec.hpp"
#include "context_modelling.hpp"
#include "mv.hpp"
#include "unit_partitioner.hpp"
#include "coding_structure.hpp"

namespace EntropyCoding {

class CABACReader {
public:
  CABACReader(BinDecoderBase &binDecoder)
      : m_BinDecoder(binDecoder), m_Bitstream(0) {}
  virtual ~CABACReader() {}

public:
  void initCtxModels(Common::Slice &slice);
  void initBitstream(Common::InputBitstream *bitstream) {
    m_Bitstream = bitstream;
    m_BinDecoder.init(m_Bitstream);
  }
  const Common::Ctx &getCtx() const { return m_BinDecoder.getCtx(); }
  Common::Ctx &getCtx() { return m_BinDecoder.getCtx(); }

public:
  // slice segment data (clause 7.3.8.1)
  bool terminating_bit();
  void remaining_bytes(bool noTrailingBytesExpected);

  // coding tree unit (clause 7.3.8.2)
  void coding_tree_unit(Common::CodingStructure &cs, const Common::UnitArea &area,
                        int (&qps)[2], unsigned ctuRsAddr);

  // sao (clause 7.3.8.3)
  void sao(Common::CodingStructure &cs, unsigned ctuRsAddr);

  void readAlfCtuFilterIndex(Common::CodingStructure &cs, unsigned ctuRsAddr);

  void ccAlfFilterControlIdc(Common::CodingStructure &cs, const Common::ComponentID compID,
                             const int curIdx, uint8_t *filterControlIdc,
                             Common::Position lumaPos, int filterCount);

  // coding (quad)tree (clause 7.3.8.4)
  void coding_tree(Common::CodingStructure &cs, Common::Partitioner &pm, Common::CUCtx &cuCtx,
                   Common::Partitioner *pPartitionerChroma = nullptr,
                   Common::CUCtx *pCuCtxChroma = nullptr);
  Common::PartSplit split_cu_mode(Common::CodingStructure &cs, Common::Partitioner &pm);
  Common::ModeType mode_constraint(Common::CodingStructure &cs, Common::Partitioner &pm,
                           const Common::PartSplit splitMode);

  // coding unit (clause 7.3.8.5)
  void coding_unit(Common::CodingUnit &cu, Common::Partitioner &pm, Common::CUCtx &cuCtx);
  void cu_skip_flag(Common::CodingUnit &cu);
  void pred_mode(Common::CodingUnit &cu);
  void bdpcm_mode(Common::CodingUnit &cu, const Common::ComponentID compID);
  void cu_pred_data(Common::CodingUnit &cu);
  void cu_bcw_flag(Common::CodingUnit &cu);
  void extend_ref_line(Common::CodingUnit &cu);
  void intra_luma_pred_modes(Common::CodingUnit &cu);
  void intra_chroma_pred_modes(Common::CodingUnit &cu);
  bool intra_chroma_lmc_mode(Common::PredictionUnit &pu);
  void intra_chroma_pred_mode(Common::PredictionUnit &pu);
  void cu_residual(Common::CodingUnit &cu, Common::Partitioner &pm, Common::CUCtx &cuCtx);
  void rqt_root_cbf(Common::CodingUnit &cu);
  void adaptive_color_transform(Common::CodingUnit &cu);
  void sbt_mode(Common::CodingUnit &cu);
  void end_of_ctu(Common::CodingUnit &cu, Common::CUCtx &cuCtx);
  void mip_flag(Common::CodingUnit &cu);
  void mip_pred_modes(Common::CodingUnit &cu);
  void mip_pred_mode(Common::PredictionUnit &pu);
  void cu_palette_info(Common::CodingUnit &cu, Common::ComponentID compBegin, uint32_t numComp,
                       Common::CUCtx &cuCtx);
  void cuPaletteSubblockInfo(Common::CodingUnit &cu, Common::ComponentID compBegin,
                             uint32_t numComp, int subSetId,
                             uint32_t &prevRunPos, unsigned &prevRunType);
  // prediction unit (clause 7.3.8.6)
  void prediction_unit(Common::PredictionUnit &pu, Common::MergeCtx &mrgCtx);
  void merge_flag(Common::PredictionUnit &pu);
  void merge_data(Common::PredictionUnit &pu);
  void affine_flag(Common::CodingUnit &cu);
  void subblock_merge_flag(Common::CodingUnit &cu);
  void merge_idx(Common::PredictionUnit &pu);
  void mmvd_merge_idx(Common::PredictionUnit &pu);
  void imv_mode(Common::CodingUnit &cu, Common::MergeCtx &mrgCtx);
  void affine_amvr_mode(Common::CodingUnit &cu, Common::MergeCtx &mrgCtx);
  void inter_pred_idc(Common::PredictionUnit &pu);
  void ref_idx(Common::PredictionUnit &pu, Common::RefPicList eRefList);
  void mvp_flag(Common::PredictionUnit &pu, Common::RefPicList eRefList);
  void Ciip_flag(Common::PredictionUnit &pu);
  void smvd_mode(Common::PredictionUnit &pu);

  // transform tree (clause 7.3.8.8)
  void transform_tree(Common::CodingStructure &cs, Common::Partitioner &pm, Common::CUCtx &cuCtx,
                      const Common::PartSplit ispType = Common::TU_NO_ISP,
                      const int subTuIdx = -1);
  bool cbf_comp(Common::CodingStructure &cs, const Common::CompArea &area, unsigned depth,
                const bool prevCbf = false, const bool useISP = false);

  // mvd coding (clause 7.3.8.9)
  void mvd_coding(Common::Mv &rMvd);

  // transform unit (clause 7.3.8.10)
  void transform_unit(Common::TransformUnit &tu, Common::CUCtx &cuCtx, Common::Partitioner &pm,
                      const int subTuCounter = -1);
  void cu_qp_delta(Common::CodingUnit &cu, int predQP, int8_t &qp);
  void cu_chroma_qp_offset(Common::CodingUnit &cu);

  // residual coding (clause 7.3.8.11)
  void residual_coding(Common::TransformUnit &tu, Common::ComponentID compID, Common::CUCtx &cuCtx);
  void ts_flag(Common::TransformUnit &tu, Common::ComponentID compID);
  void mts_idx(Common::CodingUnit &cu, Common::CUCtx &cuCtx);
  void residual_lfnst_mode(Common::CodingUnit &cu, Common::CUCtx &cuCtx);
  void isp_mode(Common::CodingUnit &cu);
  int last_sig_coeff(Common::CoeffCodingContext &cctx, Common::TransformUnit &tu,
                     Common::ComponentID compID);
  void residual_coding_subblock(Common::CoeffCodingContext &cctx, Common::TCoeff *coeff,
                                const int stateTransTable, int &state);
  void residual_codingTS(Common::TransformUnit &tu, Common::ComponentID compID);
  void residual_coding_subblockTS(Common::CoeffCodingContext &cctx, Common::TCoeff *coeff,
                                  int riceParam);
  void joint_cb_cr(Common::TransformUnit &tu, const int cbfMask);

private:
  unsigned unary_max_symbol(unsigned ctxId0, unsigned ctxIdN,
                            unsigned maxSymbol);
  unsigned unary_max_eqprob(unsigned maxSymbol);
  unsigned exp_golomb_eqprob(unsigned count);
  unsigned get_num_bits_read() { return m_BinDecoder.getNumBitsRead(); }

  void xReadTruncBinCode(uint32_t &symbol, uint32_t maxSymbol);
  void parseScanRotationModeFlag(Common::CodingUnit &cu, Common::ComponentID compBegin);
  void xDecodePLTPredIndicator(Common::CodingUnit &cu, uint32_t maxPLTSize,
                               Common::ComponentID compBegin);
  void xAdjustPLTIndex(Common::CodingUnit &cu, Common::Pel curLevel, uint32_t idx,
                       Common::PelBuf &paletteIdx, Common::PLTtypeBuf &paletteRunType,
                       int maxSymbol, Common::ComponentID compBegin);

public:
private:
  BinDecoderBase &m_BinDecoder;
  Common::InputBitstream *m_Bitstream;
  Common::ScanElement *m_scanOrder;
};

class CABACDecoder {
public:
  CABACDecoder()
      : m_CABACReaderStd(m_BinDecoderStd), m_CABACReader{&m_CABACReaderStd} {}

  CABACReader *getCABACReader(int id) { return m_CABACReader[id]; }

private:
  BinDecoder_Std m_BinDecoderStd;
  CABACReader m_CABACReaderStd;
  CABACReader *m_CABACReader[Common::BPM_NUM - 1];
};
} // namespace EntropyCoding

#endif
