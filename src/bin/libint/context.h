
#include <entity.h>

#ifndef _libint2_src_bin_libint_codecontext_h_
#define _libint2_src_bin_libint_codecontext_h_

namespace libint2 {

  /**
     CodeContext provides context for generating code
  */
  class CodeContext {
  public:
    virtual ~CodeContext() {}

    /// turn comments on and off (the default is on)
    void set_comment(bool on);
    /// returns true if to print comments
    bool comments_on() const;

    /// this function resets the context to be used for the next source file
    virtual void reset();

    /// std_header() returns declarations necessary for every source file
    virtual std::string std_header() const =0;
    /// label_to_name(label) converts label to a name valid within the context of the language
    virtual std::string label_to_name(const std::string& label) const =0;
    /// comment(statement) returns commented statement
    virtual std::string comment(const std::string& statement) const =0;
    /// open a code block
    virtual std::string open_block() const =0;
    /// close a code block
    virtual std::string close_block() const =0;

    /// unique_name<T> returns a unique name for a variable of type T
    template <typename T>
      std::string unique_name();
    /// type_name<T> returns name for type T
    template <typename T>
      std::string type_name();

  protected:
    CodeContext();
    /// generates a unique name for a floating-point variable
    virtual std::string unique_fp_name() =0;
    /// generates a unique name for an integer
    virtual std::string unique_int_name() =0;

    /// next fp index
    unsigned int next_fp_index();
    /// next int index
    unsigned int next_int_index();
    /// replaces every appearance of From with To in S
    static std::string replace_chars(const std::string& S,
                                     const std::string& From,
                                     const std::string& To);

    /// returns name of void type
    virtual std::string void_type() const =0;
    /// returns name of integer type
    virtual std::string int_type() const =0;
    /// returns name of floating-point type
    virtual std::string fp_type() const =0;
    /// returns name of pointer to floating-point type
    virtual std::string ptr_fp_type() const =0;

  private:
    unsigned int next_index_[EntityTypes::ntypes];
    void zero_out_counters();
    bool comments_on_;

  };


  /**
     CppCodeContext is an implementation of CodeContext for C++
  */
  class CppCodeContext : public CodeContext {
  public:
    CppCodeContext();
    ~CppCodeContext();

    /// Implementation of CodeContext::std_header()
    std::string std_header() const;
    /// Implementation of CodeContext::label_to_name(label)
    std::string label_to_name(const std::string& label) const;
    /// Implementation of CodeContext::comment(statement)
    std::string comment(const std::string& statement) const;
    /// Implementation of CodeContext::open_block()
    std::string open_block() const;
    /// Implementation of CodeContext::close_block()
    std::string close_block() const;

  private:
    /// Implementation of CodeContext::unique_fp_name()
    std::string unique_fp_name();
    /// Implementation of CodeContext::unique_int_name()
    std::string unique_int_name();

    std::string void_type() const;
    std::string int_type() const;
    std::string fp_type() const;
    std::string ptr_fp_type() const;

  };

};

#endif
