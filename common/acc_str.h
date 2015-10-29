//////////////////////////////////////////////////////////////////////////////
//// This software module is developed by SciDM (Scientific Data Management) in 1998-2015
//// 
//// This program is free software; you can redistribute, reuse,
//// or modify it with no restriction, under the terms of the MIT License.
//// 
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//// 
//// For any questions please contact Denis Kaznadzey at dkaznadzey@yahoo.com
//////////////////////////////////////////////////////////////////////////////

#ifndef __acc_str_h__
#define __acc_str_h__

#include <common_str.h>

#ifndef __acc_str_cpp__
#if 0
extern const char* TRUE_STR;
extern const char* FALSE_STR;
extern const char* YES_STR;
extern const char* NO_STR;
extern const char* Y_STR;
extern const char* N_STR;
extern const char* T_STR;
extern const char* F_STR;
extern const char* EMPTY_STR;
extern const char* ZERO_STR;
extern const char* ONE_STR;
extern const char* MINUS_ONE_STR;
extern const char* TWO_STR;
extern const char* MINUS_TWO_STR;
extern const char* THREE_STR;
extern const char* MINUS_THREE_STR;
extern const char* NAME_STR;
extern const char* NUMBER_STR;
extern const char* FILENAME_STR;
extern const char* INTEGER_STR;
extern const char* BOOLEAN_STR;
extern const char* STRING_STR;
extern const char* FLOAT_STR;
extern const char* DOUBLE_STR;
extern const char* OBJNAME_STR;
#endif

extern const char* FAIL_OPEN_SUBSYSTEM;
extern const char* FAIL_CLONE_TYPE;
extern const char* FAIL_NARROW_TYPE;
extern const char* FAIL_NARROW_IDSETMGR;
extern const char* FAIL_NARROW_LOCKMGR;
extern const char* FAIL_NARROW_RIGHTSMGR;
extern const char* FAIL_NARROW_SET;
extern const char* FAIL_OPEN_ACCESSOR;
extern const char* FAIL_NARROW_ACCESSOR;

extern const char* TEXT_TYPE_NAME;
extern const char* DICTIONARY_TYPE_NAME;
extern const char* ORDER_TYPE_NAME;

extern const char* OBJNUMBERING_TYPE_NAME;

extern const char* INFOSOURCE_TYPE_NAME;
extern const char* CALCLOG_TYPE_NAME;
extern const char* UPDATE_TYPE_NAME;
extern const char* ALIAS_TYPE_NAME;
extern const char* NOMENCLATURE_TYPE_NAME;
extern const char* SYNONYM_TYPE_NAME;

extern const char* CLASSIFICATION_TYPE_NAME;
extern const char* CLASSIFIEDTYPES_TYPE_NAME;
extern const char* HIERARCHYLEVEL_TYPE_NAME;
extern const char* CLASSINSTANCE_TYPE_NAME;
extern const char* ANNOTATION_TYPE_NAME;
extern const char* EVIDENCE_TYPE_NAME;
extern const char* ASSOCIATION_TYPE_NAME;
extern const char* PAIR_TYPE_NAME;
extern const char* CLUSTERING_TYPE_NAME;
extern const char* CLUSTER_TYPE_NAME;
extern const char* CLUSTER_MAP_TYPE_NAME;
extern const char* DAGVERTEX_TYPE_NAME;
extern const char* DAGEDGE_TYPE_NAME;
extern const char* MWCLUSTER_TYPE_NAME;
extern const char* MWCINDEX_TYPE_NAME;

extern const char* PHILODISTMATRIX_TYPE_NAME;
extern const char* GENETICCODE_TYPE_NAME;
extern const char* AMINOACID_TYPE_NAME;

extern const char* LITREF_TYPE_NAME;
extern const char* ORGANISM_TYPE_NAME;

extern const char* PSEQ_TYPE_NAME;
extern const char* NSEQ_TYPE_NAME;
extern const char* PSEQSET_TYPE_NAME;
extern const char* NSEQSET_TYPE_NAME;
extern const char* SEQINFO_TYPE_NAME;
extern const char* FEATURE_TYPE_NAME;
extern const char* FEATUREDATA_TYPE_NAME;
extern const char* FEATURELOCATIONMAP_TYPE_NAME;
extern const char* FRAGSET_TYPE_NAME;
extern const char* SEQLIT_TYPE_NAME;

extern const char* KTUPLECOUNTS_TYPE_NAME;
extern const char* WEIGHTMATRIX_TYPE_NAME;
extern const char* SIMMATRIX_TYPE_NAME;
extern const char* SIM_TYPE_NAME;
extern const char* SIMBATCHES_TYPE_NAME;

extern const char* SEQMODEL_TYPE_NAME;
extern const char* DOMAINLEAF_TYPE_NAME;
extern const char* DOMAINBRANCH_TYPE_NAME;

extern const char* ALIGNMENT_TYPE_NAME;
extern const char* BLOCKELEM_TYPE_NAME;
extern const char* BLOCK_TYPE_NAME;
extern const char* BLOCKALIGNMENT_TYPE_NAME;

extern const char* CLONELIBRARY_TYPE_NAME;
extern const char* CLONE_TYPE_NAME;
extern const char* CLONESEQ_TYPE_NAME;

extern const char* ACCESSION_NOMEN_ABBR;
extern const char* VERSION_NOMEN_ABBR;
extern const char* LOCUS_NOMEN_ABBR;
extern const char* GID_NOMEN_ABBR;
extern const char* NCBI_TAXID_NOMEN_ABBR;

extern const char* FEATURE_KEYS_CLASSIFICATION_NAME;
extern const char* FEATURE_QUALS_CLASSIFICATION_NAME;
extern const char* NCBI_TAXONOMY_CLASSIFICATION_NAME;
extern const char* EBI_TAXONOMY_CLASSIFICATION_NAME;
extern const char* SEQ_FILE_PATH_CLASSIFICATION_NAME;
extern const char* SEQ_FILE_PATH_CLASSIFICATION_DESCR;


extern const char* NUCL_GC_ASSOC_NAME;
extern const char* MITO_GC_ASSOC_NAME;
extern const char* CHLO_GC_ASSOC_NAME;

extern const char* SYSTEM_TYPE_CLASSIF_NAME;
extern const char* SYSTEM_TYPE_CLASSIF_DESC;

extern const char* NUCL_GC_CLASS_NAME;
extern const char* MITO_GC_CLASS_NAME;
extern const char* CHLO_GC_CLASS_NAME;

extern const char* SEQ_SIM_CLUSTERING_NAME;
extern const char* SEQ_MAP_CLUSTERING_NAME;


#endif

#endif // __seqtools_acc_str_h__
