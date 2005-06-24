
#include <ostream>
#include <string>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_defaultparams_h_
#define _libint2_src_bin_libint_defaultparams_h_

/**
  Defaults definitions for various parameters assumed by Libint
*/

namespace libint2 {
  
  class CompilationParameters {
    public:
    /// Use default parameters
    CompilationParameters();
    ~CompilationParameters();
    
    /// returns max AM for ERI
    unsigned int max_am_eri() const {
      return max_am_eri_;
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
    /// returns whether to use C-style linking
    bool use_C_linking() const {
      return use_C_linking_;
    }
    
    /// set max AM for ERI
    void max_am_eri(unsigned int a) {
      max_am_eri_ = a;
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
    /// set whether to use C style linking
    void use_C_linking(bool a) {
      use_C_linking_ = a;
    }
    
    /// print params out
    void print(std::ostream& os) const;
    
    private:
    struct Defaults {
      /// By default compile for p-functions
      static const unsigned int max_am_eri = 1;
      /// Do not vectorize by default
      static const unsigned int max_vector_length = 1;
      /// Vectorize all body by default
      static const bool vectorize_by_line = false;
      /// Produce quartet-level code by default
      static const unsigned int unroll_threshold = 1;
      /// Where to put generated library source
      static const std::string source_directory;
      /// Use C-style linking convention by default
      static const bool use_C_linking = true;
    };
    
    /// max AM for ERI
    unsigned int max_am_eri_;
    /// max vector length
    unsigned int max_vector_length_;
    /// whether to vectorize line-by-line
    bool vectorize_by_line_;
    /// unroll threshold
    unsigned int unroll_threshold_;
    /// source directory
    std::string source_directory_;
    /// whether to use C linking
    bool use_C_linking_;
  };
  
  /** The class maintains various parameters that are computed in the process of compilation
      (max stack size, etc.). Obviously, this is a Singleton.
    */
  class LibraryParameters {
    public:
    static LibraryParameters& get_library_params();
    ~LibraryParameters() {}
    
    /// returns max stack size
    unsigned int max_stack_size() const {
      return max_stack_size_;
    }
    
    /// if max_stack_size_ < size then set max_stack_size_=size
    void max_stack_size(unsigned int size) {
      if (max_stack_size_ < size)
        max_stack_size_ = size;
    }
    
    private:
    LibraryParameters();
    
    unsigned int max_stack_size_;
    
    static LibraryParameters LP_obj_;
  };
  
  struct StaticDefinitions {
    /// De facto am limit
    static const unsigned int num_am_letters = 22;
    /// am -> char conversion
    static const char am_letters[num_am_letters];
  };

  /// Converts a computation label to the name of the function
  std::string label_to_funcname(const std::string& label);

};

#endif

