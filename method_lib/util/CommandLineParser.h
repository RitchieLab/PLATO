/*
 * CommandLineParser.hpp
 *
 *  Created on: Nov 20, 2013
 *      Author: jrw32
 */

#ifndef UTILITY_COMMANDLINEPARSER_H
#define UTILITY_COMMANDLINEPARSER_H

#include <boost/program_options.hpp>

#include <boost/detail/workaround.hpp>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <vector>
#include <string>

namespace PLATO{
namespace Utility{


template<class charT>
class custom_command_line_parser : public boost::program_options::detail::cmdline {
public:
    /** Creates a command line parser for the specified arguments
        list. The 'args' parameter should not include program name.
    */
    custom_command_line_parser(const std::vector<
                                std::basic_string<charT> >& args);
    /** Creates a command line parser for the specified arguments
        list. The parameters should be the same as passed to 'main'.
    */
    custom_command_line_parser(int argc, const charT* const argv[]);

    /** Sets options descriptions to use. */
    custom_command_line_parser& options(const boost::program_options::options_description& desc);
    /** Sets positional options description to use. */
    custom_command_line_parser& positional(
        const boost::program_options::positional_options_description& desc);

    /** Sets the command line style. */
    custom_command_line_parser& style(int);
    /** Sets the extra parsers. */
    custom_command_line_parser& extra_parser(boost::program_options::ext_parser);

    /** Parses the options and returns the result of parsing.
        Throws on error.
    */
    boost::program_options::basic_parsed_options<charT> run();

    /** Specifies that unregistered options are allowed and should
        be passed though. For each command like token that looks
        like an option but does not contain a recognized name, an
        instance of basic_option<charT> will be added to result,
        with 'unrecognized' field set to 'true'. It's possible to
        collect all unrecognized options with the 'collect_unrecognized'
        funciton.
    */
    custom_command_line_parser& allow_unregistered();

    using boost::program_options::detail::cmdline::style_parser;

    custom_command_line_parser& extra_style_parser(style_parser s);

private:
    const boost::program_options::options_description* m_desc;
};

typedef custom_command_line_parser<char> CommandLineParser;
typedef custom_command_line_parser<wchar_t> UnicodeCommandLineParser;

template<class charT>
custom_command_line_parser<charT>::
custom_command_line_parser(const std::vector<
                            std::basic_string<charT> >& xargs)
    : boost::program_options::detail::cmdline(boost::program_options::to_internal(xargs))
{}


template<class charT>
custom_command_line_parser<charT>::
custom_command_line_parser(int argc, const charT* const argv[])
: boost::program_options::detail::cmdline(
    // Explicit template arguments are required by gcc 3.3.1
    // (at least mingw version), and do no harm on other compilers.
    boost::program_options::to_internal(boost::program_options::detail::make_vector<charT, const charT* const*>(argv+1, argv+argc+!argc)))
{}


template<class charT>
custom_command_line_parser<charT>&
custom_command_line_parser<charT>::options(const boost::program_options::options_description& desc)
{
    boost::program_options::detail::cmdline::set_options_description(desc);
    m_desc = &desc;
    return *this;
}

template<class charT>
custom_command_line_parser<charT>&
custom_command_line_parser<charT>::positional(
    const boost::program_options::positional_options_description& desc)
{
    boost::program_options::detail::cmdline::set_positional_options(desc);
    return *this;
}

template<class charT>
custom_command_line_parser<charT>&
custom_command_line_parser<charT>::style(int xstyle)
{
    boost::program_options::detail::cmdline::style(xstyle);
    return *this;
}

template<class charT>
custom_command_line_parser<charT>&
custom_command_line_parser<charT>::extra_parser(boost::program_options::ext_parser ext)
{
    boost::program_options::detail::cmdline::set_additional_parser(ext);
    return *this;
}

template<class charT>
custom_command_line_parser<charT>&
custom_command_line_parser<charT>::allow_unregistered()
{
    boost::program_options::detail::cmdline::allow_unregistered();
    return *this;
}

template<class charT>
custom_command_line_parser<charT>&
custom_command_line_parser<charT>::extra_style_parser(style_parser s)
{
    boost::program_options::detail::cmdline::extra_style_parser(s);
    return *this;
}

template<class charT>
boost::program_options::basic_parsed_options<charT>
custom_command_line_parser<charT>::run()
{
	using namespace boost::program_options;
	using namespace boost::program_options::command_line_style;
	using std::vector;
	using std::string;
//	using namespace boost::program_options::detail::

    // save the canonical prefixes which were used by this cmdline parser
    //    eventually inside the parsed results
    //    This will be handy to format recognisable options
    //    for diagnostic messages if everything blows up much later on
    parsed_options super_result(m_desc, detail::cmdline::get_canonical_option_prefix());

    // The parsing is done by having a set of 'style parsers'
    // and trying then in order. Each parser is passed a vector
    // of unparsed tokens and can consume some of them (by
    // removing elements on front) and return a vector of options.
    //
    // We try each style parser in turn, untill some input
    // is consumed. The returned vector of option may contain the
    // result of just syntactic parsing of token, say --foo will
    // be parsed as option with name 'foo', and the style parser
    // is not required to care if that option is defined, and how
    // many tokens the value may take.
    // So, after vector is returned, we validate them.
    assert(m_desc);

    vector<style_parser> style_parsers;

    if (m_style_parser)
        style_parsers.push_back(m_style_parser);

    if (m_additional_parser)
        style_parsers.push_back(
            boost::bind(&cmdline::handle_additional_parser, this, _1));

    if (m_style & allow_long)
        style_parsers.push_back(
            boost::bind(&cmdline::parse_long_option, this, _1));

    if ((m_style & allow_long_disguise))
        style_parsers.push_back(
            boost::bind(&cmdline::parse_disguised_long_option, this, _1));

    if ((m_style & allow_short) && (m_style & allow_dash_for_short))
        style_parsers.push_back(
            boost::bind(&cmdline::parse_short_option, this, _1));

    if ((m_style & allow_short) && (m_style & allow_slash_for_short))
        style_parsers.push_back(boost::bind(&cmdline::parse_dos_option, this, _1));

    style_parsers.push_back(boost::bind(&cmdline::parse_terminator, this, _1));

    vector<option> result;
    bool ok=true;
    while(!args.empty() && ok)
    {
        ok = false;
        for(unsigned i = 0; i < style_parsers.size(); ++i)
        {
            unsigned current_size = static_cast<unsigned>(args.size());
            vector<option> next = style_parsers[i](args);

            // Check that option names
            // are valid, and that all values are in place.
            if (!next.empty())
            {
                vector<string> e;
                for(unsigned k = 0; k < next.size()-1; ++k) {
                    finish_option(next[k], e, style_parsers);
                }
                // For the last option, pass the unparsed tokens
                // so that they can be added to next.back()'s values
                // if appropriate.
                finish_option(next.back(), args, style_parsers);
                for (unsigned j = 0; j < next.size(); ++j)
                    result.push_back(next[j]);
            }

            if (args.size() != current_size) {
                ok = true;
                break;
            }
        }

        if (!ok) {
        	while(!args.empty()){
	            option opt;
    	        opt.value.push_back(args[0]);
    	        opt.original_tokens.push_back(args[0]);
    	        result.push_back(opt);
    	        args.erase(args.begin());
    	    }
        }
    }

    /* If an key option is followed by a positional option,
        can can consume more tokens (e.g. it's multitoken option),
        give those tokens to it.  */
    vector<option> result2;
    for (unsigned i = 0; i < result.size(); ++i)
    {
        result2.push_back(result[i]);
        option& opt = result2.back();

        if (opt.string_key.empty())
            continue;

        const option_description* xd;
        try
        {
            xd = m_desc->find_nothrow(opt.string_key,
                                        is_style_active(allow_guessing),
                                        is_style_active(long_case_insensitive),
                                        is_style_active(short_case_insensitive));
        }
        catch(error_with_option_name& e)
        {
            // add context and rethrow
            e.add_context(opt.string_key, opt.original_tokens[0], get_canonical_option_prefix());
            throw;
        }

        if (!xd)
            continue;

        unsigned min_tokens = xd->semantic()->min_tokens();
        unsigned max_tokens = xd->semantic()->max_tokens();
        if (min_tokens < max_tokens && opt.value.size() < max_tokens)
        {
            // This option may grab some more tokens.
            // We only allow to grab tokens that are not already
            // recognized as key options.

            int can_take_more = max_tokens - static_cast<int>(opt.value.size());
            unsigned j = i+1;
            for (; can_take_more && j < result.size(); --can_take_more, ++j)
            {
                option& opt2 = result[j];
                if (!opt2.string_key.empty())
                    break;

                if (opt2.position_key == INT_MAX)
                {
                    // We use INT_MAX to mark positional options that
                    // were found after the '--' terminator and therefore
                    // should stay positional forever.
                    break;
                }

                assert(opt2.value.size() == 1);

                opt.value.push_back(opt2.value[0]);

                assert(opt2.original_tokens.size() == 1);

                opt.original_tokens.push_back(opt2.original_tokens[0]);
            }
            i = j-1;
        }
    }
    result.swap(result2);


    // Assign position keys to positional options.
    int position_key = 0;
    for(unsigned i = 0; i < result.size(); ++i) {
        if (result[i].string_key.empty())
            result[i].position_key = position_key++;
    }

    if (m_positional)
    {
        unsigned position = 0;
        for (unsigned i = 0; i < result.size(); ++i) {
            option& opt = result[i];
            if (opt.position_key != -1) {
                if (position >= m_positional->max_total_count())
                {
                    boost::throw_exception(too_many_positional_options_error());
                }
                opt.string_key = m_positional->name_for_position(position);
                ++position;
            }
        }
    }

    // set case sensitive flag
    for (unsigned i = 0; i < result.size(); ++i) {
        if (result[i].string_key.size() > 2 ||
                    (result[i].string_key.size() > 1 && result[i].string_key[0] != '-'))
        {
            // it is a long option
            result[i].case_insensitive = is_style_active(long_case_insensitive);
        }
        else
        {
            // it is a short option
            result[i].case_insensitive = is_style_active(short_case_insensitive);
        }
    }

    super_result.options = result;

    // Presense of parsed_options -> wparsed_options conversion
    // does the trick.
    return basic_parsed_options<charT>(super_result);
}

}
}

#endif /* COMMANDLINEPARSER_HPP_ */

