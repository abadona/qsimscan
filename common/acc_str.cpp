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

#define __acc_str_cpp__
#include "acc_str.h"

#if 0
const char* TRUE_STR = "TRUE";
const char* FALSE_STR = "FALSE";
const char* YES_STR = "YES";
const char* NO_STR = "NO";
const char* Y_STR = "Y";
const char* N_STR = "N";
const char* T_STR = "T";
const char* F_STR = "F";
const char* EMPTY_STR = "";
const char* ZERO_STR = "0";
const char* ONE_STR = "1";
const char* MINUS_ONE_STR = "-1";
const char* TWO_STR = "2";
const char* MINUS_TWO_STR = "-2";
const char* THREE_STR = "3";
const char* MINUS_THREE_STR = "-3";
const char* NAME_STR = "name";
const char* NUMBER_STR = "number";
const char* FILENAME_STR = "filename";
const char* INTEGER_STR = "integer";
const char* BOOLEAN_STR = "boolean";
const char* STRING_STR = "boolean";
const char* FLOAT_STR = "float_number";
const char* DOUBLE_STR = "double_prec_number";
const char* OBJNAME_STR = "object_name";
#endif


const char* FAIL_OPEN_SUBSYSTEM = "Subsystem open failed on Session object";
const char* FAIL_CLONE_TYPE = "clone operation failed for Type accessor";
const char* FAIL_NARROW_TYPE = "Type accessor narrowing failed";
const char* FAIL_NARROW_IDSETMGR = "IdSetMgr accessor narrowing failed";
const char* FAIL_NARROW_LOCKMGR = "Lock accessor narrowing failed";
const char* FAIL_NARROW_RIGHTSMGR = "Rights accessor narrowing failed";
const char* FAIL_NARROW_SET = "clone operation failed for Set accessor";
const char* FAIL_OPEN_ACCESSOR = "Type->openAccessor operation failed";
const char* FAIL_NARROW_ACCESSOR = "narrowing operation for Accessor failed";

const char* TEXT_TYPE_NAME = "Txt";
const char* DICTIONARY_TYPE_NAME = "Dictionary";
const char* ORDER_TYPE_NAME = "Order";

const char* OBJNUMBERING_TYPE_NAME = "ObjNumbering";

const char* INFOSOURCE_TYPE_NAME = "InfoSource";
const char* CALCLOG_TYPE_NAME = "CalcLog";
const char* UPDATE_TYPE_NAME = "Update";
const char* ALIAS_TYPE_NAME = "Alias";
const char* NOMENCLATURE_TYPE_NAME = "Nomenclature";
const char* SYNONYM_TYPE_NAME = "Synonym";

const char* CLASSIFICATION_TYPE_NAME = "Classification";
const char* CLASSIFIEDTYPES_TYPE_NAME = "ClassifiedTypes";
const char* HIERARCHYLEVEL_TYPE_NAME = "HierarchyLevel";
const char* CLASSINSTANCE_TYPE_NAME = "ClassInstance";
const char* ANNOTATION_TYPE_NAME = "Annotation";
const char* EVIDENCE_TYPE_NAME = "Evidence";
const char* ASSOCIATION_TYPE_NAME = "Association";
const char* PAIR_TYPE_NAME = "Pair";
const char* CLUSTERING_TYPE_NAME = "Clustering";
const char* CLUSTER_TYPE_NAME = "Cluster";
const char* CLUSTERMAP_TYPE_NAME = "ClusterMap";
const char* DAGVERTEX_TYPE_NAME = "DAGvertex";
const char* DAGEDGE_TYPE_NAME = "DAGedge";
const char* MWCLUSTER_TYPE_NAME = "MWCluster";
const char* MWCINDEX_TYPE_NAME = "MWCindex";

const char* PHILODISTMATRIX_TYPE_NAME = "PhiloDistMatrix";
const char* GENETICCODE_TYPE_NAME = "GeneticCode";
const char* AMINOACID_TYPE_NAME = "Aminoacid";

const char* LITREF_TYPE_NAME = "LitRef";
const char* ORGANISM_TYPE_NAME = "Organism";

const char* PSEQ_TYPE_NAME = "PSeq";
const char* NSEQ_TYPE_NAME = "NSeq";
const char* PSEQSET_TYPE_NAME = "PSeq{1}";
const char* NSEQSET_TYPE_NAME = "NSeq{1}";
const char* SEQINFO_TYPE_NAME = "SeqInfo";
const char* FEATURE_TYPE_NAME = "Feature";
const char* FEATUREDATA_TYPE_NAME = "FeatureData";
const char* FEATURELOCATIONMAP_TYPE_NAME = "FeatureLocationMap";
const char* FRAGSET_TYPE_NAME = "Fragset";
const char* SEQLIT_TYPE_NAME = "SeqLit";

const char* KTUPLECOUNTS_TYPE_NAME = "KtupleCounts";
const char* WEIGHTMATRIX_TYPE_NAME = "WeightMatrix";
const char* SIMMATRIX_TYPE_NAME = "SimMatrix";
const char* SIM_TYPE_NAME = "Sim";
const char* SIMBATCHES_TYPE_NAME = "SimBatches";

const char* SEQMODEL_TYPE_NAME = "SeqModel";
const char* DOMAINLEAF_TYPE_NAME = "DomainLeaf";
const char* DOMAINBRANCH_TYPE_NAME = "DomainBranch";

const char* ALIGNMENT_TYPE_NAME = "Alignment";
const char* BLOCKELEM_TYPE_NAME = "BlockElem";
const char* BLOCK_TYPE_NAME = "Block";
const char* BLOCKALIGNMENT_TYPE_NAME = "BlockAlignment";

const char* CLONELIBRARY_TYPE_NAME = "CloneLibrary";
const char* CLONE_TYPE_NAME = "Clone";
const char* CLONESEQ_TYPE_NAME = "CloneSeq";

const char* ACCESSION_NOMEN_ABBR = "ACC";
const char* VERSION_NOMEN_ABBR = "ACCV";
const char* LOCUS_NOMEN_ABBR = "LOC";
const char* GID_NOMEN_ABBR = "GID";
const char* NCBI_TAXID_NOMEN_ABBR = "taxid";

const char* FEATURE_KEYS_CLASSIFICATION_NAME = "FeatureKeys";
const char* FEATURE_QUALS_CLASSIFICATION_NAME = "FeatureQualifiers";
const char* NCBI_TAXONOMY_CLASSIFICATION_NAME = "NCBI_taxonomy";
const char* EBI_TAXONOMY_CLASSIFICATION_NAME = "EBI_taxonomy";
const char* SEQ_FILE_PATH_CLASSIFICATION_NAME = "SEQ_FILE_PATH";
const char* SEQ_FILE_PATH_CLASSIFICATION_DESCR = "Location of source sequence data file in the filesystem";

const char* NUCL_GC_ASSOC_NAME = "tax_nuclear_gencode";
const char* MITO_GC_ASSOC_NAME = "tax_mitochondrial_gencode";
const char* CHLO_GC_ASSOC_NAME = "tax_chloroplast_gencode";

const char* SYSTEM_TYPE_CLASSIF_NAME = "System_Type";
const char* SYSTEM_TYPE_CLASSIF_DESC = "Classification for global constant object classes";

const char* NUCL_GC_CLASS_NAME = NUCL_GC_ASSOC_NAME;
const char* MITO_GC_CLASS_NAME = MITO_GC_ASSOC_NAME;
const char* CHLO_GC_CLASS_NAME = CHLO_GC_ASSOC_NAME;

const char* SEQ_SIM_CLUSTERING_NAME = "SEQ_SIM_CLUSTERING";
const char* SEQ_MAP_CLUSTERING_NAME = "SEQ_HITMAP_CLUSTERING";

