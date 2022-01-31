#ifndef ENTROPY_CODEC_CABAC_WRITER
#define ENTROPY_CODEC_CABAC_WRITER

#include "arith_codec.hpp"
#include "bit_stream.hpp"
#include "coding_structure.hpp"
#include "context_modelling.hpp"
#include "mv.hpp"
#include "slice.hpp"
#include "unit_partitioner.hpp"

namespace EntropyCoding {

class CABACWriter {
public:
  CABACWriter(BinEncIf &binEncoder)
      : m_BinEncoder(binEncoder), m_Bitstream(0),
        m_TestCtx(m_BinEncoder.getCtx()) {}
  virtual ~CABACWriter() {}

public:
  void initCtxModels(const Common::Slice &slice);
  Common::SliceType getCtxInitId(const Common::Slice &slice);
  void initBitstream(Common::OutputBitstream *bitstream) {
    m_Bitstream = bitstream;
    m_BinEncoder.init(m_Bitstream);
  }

  const Common::Ctx &getCtx() const { return m_BinEncoder.getCtx(); }
  Common::Ctx &getCtx() { return m_BinEncoder.getCtx(); }

  void start() { m_BinEncoder.start(); }
  void resetBits() { m_BinEncoder.resetBits(); }
  uint64_t getEstFracBits() const { return m_BinEncoder.getEstFracBits(); }
  uint32_t getNumBins() { return m_BinEncoder.getNumBins(); }
  bool isEncoding() { return m_BinEncoder.isEncoding(); }

public:
  // slice segment data (clause 7.3.8.1)
  void end_of_slice();

  // coding tree unit (clause 7.3.8.2)
  void coding_tree_unit(Common::CodingStructure &cs, const Common::UnitArea &area,
                        int (&qps)[2], unsigned ctuRsAddr, bool skipSao = false,
                        bool skipAlf = false);

  // sao (clause 7.3.8.3)
  void sao(const Common::Slice &slice, unsigned ctuRsAddr);
  void sao_block_pars(const Common::SAOBlkParam &saoPars, const Common::BitDepths &bitDepths,
                      bool *sliceEnabled, bool leftMergeAvail,
                      bool aboveMergeAvail, bool onlyEstMergeInfo);
  void sao_offset_pars(const Common::SAOOffset &ctbPars, Common::ComponentID compID,
                       bool sliceEnabled, int bitDepth);
  // coding (quad)tree (clause 7.3.8.4)
  void coding_tree(const Common::CodingStructure &cs, Common::Partitioner &pm, Common::CUCtx &cuCtx,
                   Common::Partitioner *pPartitionerChroma = nullptr,
                   Common::CUCtx *pCuCtxChroma = nullptr);
  void split_cu_mode(const Common::PartSplit split, const Common::CodingStructure &cs,
                     Common::Partitioner &pm);
  void mode_constraint(const Common::PartSplit split, const Common::CodingStructure &cs,
                       Common::Partitioner &pm, const Common::ModeType modeType);

  // coding unit (clause 7.3.8.5)
  void coding_unit(const Common::CodingUnit &cu, Common::Partitioner &pm, Common::CUCtx &cuCtx);
  void cu_skip_flag(const Common::CodingUnit &cu);
  void pred_mode(const Common::CodingUnit &cu);
  void bdpcm_mode(const Common::CodingUnit &cu, const Common::ComponentID compID);

  void cu_pred_data(const Common::CodingUnit &cu);
  void cu_bcw_flag(const Common::CodingUnit &cu);
  void extend_ref_line(const Common::PredictionUnit &pu);
  void extend_ref_line(const Common::CodingUnit &cu);
  void intra_luma_pred_modes(const Common::CodingUnit &cu);
  void intra_luma_pred_mode(const Common::PredictionUnit &pu);
  void intra_chroma_pred_modes(const Common::CodingUnit &cu);
  void intra_chroma_lmc_mode(const Common::PredictionUnit &pu);
  void intra_chroma_pred_mode(const Common::PredictionUnit &pu);
  void cu_residual(const Common::CodingUnit &cu, Common::Partitioner &pm, Common::CUCtx &cuCtx);
  void rqt_root_cbf(const Common::CodingUnit &cu);
  void adaptive_color_transform(const Common::CodingUnit &cu);
  void sbt_mode(const Common::CodingUnit &cu);
  void end_of_ctu(const Common::CodingUnit &cu, Common::CUCtx &cuCtx);
  void mip_flag(const Common::CodingUnit &cu);
  void mip_pred_modes(const Common::CodingUnit &cu);
  void mip_pred_mode(const Common::PredictionUnit &pu);
  void cu_palette_info(const Common::CodingUnit &cu, Common::ComponentID compBegin,
                       uint32_t numComp, Common::CUCtx &cuCtx);
  void cuPaletteSubblockInfo(const Common::CodingUnit &cu, Common::ComponentID compBegin,
                             uint32_t numComp, int subSetId,
                             uint32_t &prevRunPos, unsigned &prevRunType);
  Common::Pel writePLTIndex(const Common::CodingUnit &cu, uint32_t idx, Common::PelBuf &paletteIdx,
                    Common::PLTtypeBuf &paletteRunType, int maxSymbol,
                    Common::ComponentID compBegin);
  // prediction unit (clause 7.3.8.6)
  void prediction_unit(const Common::PredictionUnit &pu);
  void merge_flag(const Common::PredictionUnit &pu);
  void merge_data(const Common::PredictionUnit &pu);
  void affine_flag(const Common::CodingUnit &cu);
  void subblock_merge_flag(const Common::CodingUnit &cu);
  void merge_idx(const Common::PredictionUnit &pu);
  void mmvd_merge_idx(const Common::PredictionUnit &pu);
  void imv_mode(const Common::CodingUnit &cu);
  void affine_amvr_mode(const Common::CodingUnit &cu);
  void inter_pred_idc(const Common::PredictionUnit &pu);
  void ref_idx(const Common::PredictionUnit &pu, Common::RefPicList eRefList);
  void mvp_flag(const Common::PredictionUnit &pu, Common::RefPicList eRefList);

  void Ciip_flag(const Common::PredictionUnit &pu);
  void smvd_mode(const Common::PredictionUnit &pu);

  // transform tree (clause 7.3.8.8)
  void transform_tree(const Common::CodingStructure &cs, Common::Partitioner &pm, Common::CUCtx &cuCtx,
                      const Common::PartSplit ispType = Common::TU_NO_ISP,
                      const int subTuIdx = -1);
  void cbf_comp(const Common::CodingStructure &cs, bool cbf, const Common::CompArea &area,
                unsigned depth, const bool prevCbf = false,
                const bool useISP = false);

  // mvd coding (clause 7.3.8.9)
  void mvd_coding(const Common::Mv &rMvd, int8_t imv);
  // transform unit (clause 7.3.8.10)
  void transform_unit(const Common::TransformUnit &tu, Common::CUCtx &cuCtx, Common::Partitioner &pm,
                      const int subTuCounter = -1);
  void cu_qp_delta(const Common::CodingUnit &cu, int predQP, const int8_t qp);
  void cu_chroma_qp_offset(const Common::CodingUnit &cu);

  // residual coding (clause 7.3.8.11)
  void residual_coding(const Common::TransformUnit &tu, Common::ComponentID compID,
                       Common::CUCtx *cuCtx = nullptr);
  void ts_flag(const Common::TransformUnit &tu, Common::ComponentID compID);
  void mts_idx(const Common::CodingUnit &cu, Common::CUCtx *cuCtx);
  void residual_lfnst_mode(const Common::CodingUnit &cu, Common::CUCtx &cuCtx);
  void isp_mode(const Common::CodingUnit &cu);
  void last_sig_coeff(Common::CoeffCodingContext &cctx, const Common::TransformUnit &tu,
                      Common::ComponentID compID);
  void residual_coding_subblock(Common::CoeffCodingContext &cctx, const Common::TCoeff *coeff,
                                const int stateTransTable, int &state);
  void residual_codingTS(const Common::TransformUnit &tu, Common::ComponentID compID);
  void residual_coding_subblockTS(Common::CoeffCodingContext &cctx, const Common::TCoeff *coeff,
                                  unsigned (&RiceBit)[8], int riceParam,
                                  bool ricePresentFlag);
  void joint_cb_cr(const Common::TransformUnit &tu, const int cbfMask);

  void codeAlfCtuEnableFlags(const Common::CodingStructure &cs, Common::ChannelType channel,
                             Common::AlfParam *alfParam);
  void codeAlfCtuEnableFlags(const Common::CodingStructure &cs, Common::ComponentID compID,
                             Common::AlfParam *alfParam);
  void codeAlfCtuEnableFlag(const Common::CodingStructure &cs, uint32_t ctuRsAddr,
                            const int compIdx, Common::AlfParam *alfParam);
  void codeAlfCtuFilterIndex(const Common::CodingStructure &cs, uint32_t ctuRsAddr,
                             bool alfEnableLuma);

  void codeAlfCtuAlternatives(const Common::CodingStructure &cs, Common::ChannelType channel,
                              Common::AlfParam *alfParam);
  void codeAlfCtuAlternatives(const Common::CodingStructure &cs, Common::ComponentID compID,
                              Common::AlfParam *alfParam);
  void codeAlfCtuAlternative(const Common::CodingStructure &cs, uint32_t ctuRsAddr,
                             const int compIdx,
                             const Common::AlfParam *alfParam = NULL);
  void codeCcAlfFilterControlIdc(uint8_t idcVal, const Common::CodingStructure &cs,
                                 const Common::ComponentID compID, const int curIdx,
                                 const uint8_t *filterControlIdc,
                                 Common::Position lumaPos, const int filterCount);

private:
  void unary_max_symbol(unsigned symbol, unsigned ctxId0, unsigned ctxIdN,
                        unsigned maxSymbol);
  void unary_max_eqprob(unsigned symbol, unsigned maxSymbol);
  void exp_golomb_eqprob(unsigned symbol, unsigned count);

  // statistic
  unsigned get_num_written_bits() { return m_BinEncoder.getNumWrittenBits(); }

  void xWriteTruncBinCode(uint32_t uiSymbol, uint32_t uiMaxSymbol);
  void codeScanRotationModeFlag(const Common::CodingUnit &cu, Common::ComponentID compBegin);
  void xEncodePLTPredIndicator(const Common::CodingUnit &cu, uint32_t maxPltSize,
                               Common::ComponentID compBegin);

private:
  BinEncIf &m_BinEncoder;
  Common::OutputBitstream *m_Bitstream;
  Common::Ctx m_TestCtx;
  Common::ScanElement *m_scanOrder;
};

class CABACEncoder {
public:
  CABACEncoder()
      : m_CABACWriterStd(m_BinEncoderStd),
        m_CABACEstimatorStd(m_BitEstimatorStd),
        m_CABACWriter{
            &m_CABACWriterStd,
        },
        m_CABACEstimator{&m_CABACEstimatorStd} {}

  CABACWriter *getCABACWriter(const Common::SPS *sps) { return m_CABACWriter[0]; }
  CABACWriter *getCABACEstimator(const Common::SPS *sps) { return m_CABACEstimator[0]; }

private:
  BinEncoder_Std m_BinEncoderStd;
  BitEstimator_Std m_BitEstimatorStd;
  CABACWriter m_CABACWriterStd;
  CABACWriter m_CABACEstimatorStd;
  CABACWriter *m_CABACWriter[Common::BPM_NUM - 1];
  CABACWriter *m_CABACEstimator[Common::BPM_NUM - 1];
};
} // namespace EntropyCoding

#endif //ENTROPY_CODEC_CABAC_WRITER
