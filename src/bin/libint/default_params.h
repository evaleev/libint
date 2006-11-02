
#include <ostream>
#include <string>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_defaultparams_h_
#define _libint2_src_bin_libint_defaultparams_h_

/**
  Defaults definitions for various parameters assumed by Libint
*/

namespace libint2 {

  /// These are the parameters received by the compiler
  class CompilationParameters {
    public:
    /// Use default parameters
    CompilationParameters();
    ~CompilationParameters();
    
    /// returns max AM for general integrals
    unsigned int max_am() const {
      return max_am_;
    }
    /// returns max AM of general integrals for which to produce optimal code
    unsigned int max_am_opt() const {
      return max_am_opt_;
    }
    /// returns max AM for ERI
    unsigned int max_am_eri() const {
      return max_am_eri_;
    }
    /// returns max AM of ERI for which to produce optimal code
    unsigned int max_am_eri_opt() const {
      return max_am_eri_opt_;
    }
    /// returns max AM for G12 integrals
    unsigned int max_am_g12() const {
      return max_am_g12_;
    }
    /// returns max AM of G12 ints for which to produce optimal code
    unsigned int max_am_g12_opt() const {
      return max_am_g12_opt_;
    }
    /// returns max vector length
    unsigned int max_vector_length() const {
      return max_vector_length_;
    }
    /// returns whether to vectorize line-by-line
    bool vectorize_by_line() const {
      return vectorize_by_line_;
    }
    /// returns unroll threshold
    unsigned int unroll_threshold() const {
      return unroll_threshold_;
    }
    /// returns directory path for the generated source
    const std::string& source_directory() const {
      return source_directory_;
    }
    /// the API prefix
    const std::string& api_prefix() const {
      return api_prefix_;
    }
    /// returns whether to use C-style linking
    bool use_C_linking() const {
      return use_C_linking_;
    }
    /// whether FLOP counting is enabled
    bool count_flops() const {
      return count_flops_;
    }
    /// whether target integrals are accumulated
    bool accumulate_targets() const {
      return accumulate_targets_;
    }
    /// name of the floating-point type
    const std::string& realtype() const {
      return realtype_;
    }
    
    /// set max AM for general integrals
    void max_am(unsigned int a) {
      max_am_ = a;
    }
    /// set max AM for "optimized" integrals
    void max_am_opt(unsigned int a) {
      max_am_opt_ = a;
    }
    /// set max AM for ERI
    void max_am_eri(unsigned int a) {
      max_am_eri_ = a;
    }
    /// set max AM for "optimized" ERI
    void max_am_eri_opt(unsigned int a) {
      max_am_eri_opt_ = a;
    }
    /// set max AM for G12
    void max_am_g12(unsigned int a) {
      max_am_g12_ = a;
    }
    /// set max AM for "optimized" G12
    void max_am_g12_opt(unsigned int a) {
      max_am_g12_opt_ = a;
    }
    /// set max vector length
    void max_vector_length(unsigned int a) {
      max_vector_length_ = a;
    }
    /// set vectorize_by_line flag
    void vectorize_by_line(bool flag) {
      vectorize_by_line_ = flag;
    }
    /// set unroll threshold
    void unroll_threshold(unsigned int a) {
      unroll_threshold_ = a;
    }
    /// set generated source directory
    void source_directory(const std::string& a) {
      source_directory_ = a;
    }
    /// API prefix
    void api_prefix(const std::string& a) {
      api_prefix_ = a;
    }
    /// set whether to use C style linking
    void use_C_linking(bool a) {
      use_C_linking_ = a;
    }
    /// set whether to count FLOPs
    void count_flops(bool a) {
      count_flops_ = a;
    }
    /// accumulate targets?
    void accumulate_targets(bool a) {
      accumulate_targets_ = a;
    }
    /// which floating-point type to use
    void realtype(const std::string& realtype) {
      realtype_ = realtype;
    }
    
    /// print params out
    void print(std::ostream& os) const;
    
    private:
    struct Defaults {
      /// By default compile general integrals for p-functions
      static const unsigned int max_am = 1;
      /// By default compile ERI for p-functions
      static const unsigned int max_am_eri = 1;
      /// By default compile G12 integrals for p-functions
      static const unsigned int max_am_g12 = 1;
      /// Do not vectorize by default
      static const unsigned int max_vector_length = 1;
      /// Vectorize all body by default
      static const bool vectorize_by_line = false;
      /// Produce quartet-level code by default
      static const unsigned int unroll_threshold = 1;
      /// Where to put generated library source
      static const std::string source_directory;
      /// API prefix
      static const std::string api_prefix;
      /// Use C-style linking convention by default
      static const bool use_C_linking = true;
      /// Do not count FLOPs by default
      static const bool count_flops = false;
      /// Do not accumulate targets by default
      static const bool accumulate_targets = false;
      /// Use double for computations
      static const std::string realtype;
    };
    
    /// max AM for general integrals
    unsigned int max_am_;
    /// max AM for "optimized" general integrals
    unsigned int max_am_opt_;
    /// max AM for ERI
    unsigned int max_am_eri_;
    /// max AM for "optimized" ERI
    unsigned int max_am_eri_opt_;
    /// max AM for G12
    unsigned int max_am_g12_;
    /// max AM for "optimized" G12
    unsigned int max_am_g12_opt_;
    /// max vector length
    unsigned int max_vector_length_;
    /// whether to vectorize line-by-line
    bool vectorize_by_line_;
    /// unroll threshold
    unsigned int unroll_threshold_;
    /// source directory
    std::string source_directory_;
    /// API prefix
    std::string api_prefix_;
    /// whether to use C linking
    bool use_C_linking_;
    /// whether to count FLOPs
    bool count_flops_;
    /// whether to accumulate targets
    bool accumulate_targets_;
    /// name of the floating-point type
    std::string realtype_;
  };
  
  /** This class maintains various parameters for each task type
      which can only be determined during the source generation
      (max stack size, etc.).
    */
  class TaskParameters {
    public:
    TaskParameters();
    ~TaskParameters() {}
    
    /// returns the max number of targets
    unsigned int max_ntarget() const {
      return max_ntarget_;
    }
    /// returns max stack size
    unsigned int max_stack_size() const {
      return max_stack_size_;
    }
    /** returns max vector stack size.
        vector stack is only used to hold intermediate quantities
        in set-level RR code. This is only needed when doing linewise vectorization.
      */
    unsigned int max_vector_stack_size() const {
      return max_vector_stack_size_;
    }
    /** returns max rank of high-significance functions in a HRR call.
        This is only needed when doing linewise vectorization.
      */
    unsigned int max_hrr_hsrank() const {
      return max_hrr_hsrank_;
    }
    /** returns max rank of low-significance functions in a HRR call.
        This is only needed when doing linewise vectorization.
      */
    unsigned int max_hrr_lsrank() const {
      return max_hrr_lsrank_;
    }
    
    /// if max_ntarget_ < ntarget then set max_ntarget_=ntarget
    void max_ntarget(unsigned int ntarget) {
      if (max_ntarget_ < ntarget)
        max_ntarget_ = ntarget;
    }
    
    /// if max_stack_size_ < size then set max_stack_size_=size
    void max_stack_size(unsigned int size) {
      if (max_stack_size_ < size)
        max_stack_size_ = size;
    }
    
    /// if max_vector_stack_size_ < size then set max_vector_stack_size_=size
    void max_vector_stack_size(unsigned int size) {
      if (max_vector_stack_size_ < size)
        max_vector_stack_size_ = size;
    }

    /// if max_hrr_hsrank_ < rank then set max_hrr_hsrank_=rank
    void max_hrr_hsrank(unsigned int rank) {
      if (max_hrr_hsrank_ < rank)
        max_hrr_hsrank_ = rank;
    }
    
    /// if max_hrr_lsrank_ < rank then set max_hrr_lsrank_=rank
    void max_hrr_lsrank(unsigned int rank) {
      if (max_hrr_lsrank_ < rank)
        max_hrr_lsrank_ = rank;
    }
    
    private:
    unsigned int max_ntarget_;
    unsigned int max_stack_size_;
    unsigned int max_vector_stack_size_;
    unsigned int max_hrr_hsrank_;
    unsigned int max_hrr_lsrank_;
  };

  /// Static parameters
  struct StaticDefinitions {
    /// De facto am limit
    static const unsigned int num_am_letters = 22;
    /// am -> char conversion
    static const char am_letters[num_am_letters];
  };

  /// Converts a label, e.g. name of the target node, to the name of the function to compute it
  std::string label_to_funcname(const std::string& label);

};

#endif

