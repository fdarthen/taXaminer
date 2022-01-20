"""Modify HTML

Modifies the HTML 3D plot file to make is self-contained and
the text in the hoverwindow selectable.

Expects prepared config file to find path to HTML file
(preparation of config by prepare_and_check.py)
"""

import logging
from bs4 import BeautifulSoup as bs
from jsmin import jsmin
import sys


def make_selectable(script_path):
    """Make text in hoverwindows in 3D plot selectable.

    Replaces parts of javascript code to make the text in the
    hoverwindow selectable. Highly specified function and may not work
    with other versions of htmlwidgets and plotly (developed for
    htmlwidgets-1.5.7 and plotly-4.9.0; compatiblity tested with
    htmlwidgets-1.5.4 and plotly-4.10.0)

    Args:
      script_path(str): path to javascript file where to replace code
    """

    with open(script_path, 'r') as file_in:
        filedata = file_in.read()

    old = ('r._paperdiv=r._container.selectAll(".svg-container").data([0]),'
        'r._paperdiv.enter().append("div").classed("user-select-none",!0).classed("svg-container",!0).style("position","relative"),'
        'r._glcontainer=r._paperdiv.selectAll(".gl-container").data([{}]),r._glcontainer.enter().append("div").classed("gl-container",!0),'
        'r._paperdiv.selectAll(".main-svg").remove(),r._paperdiv.select(".modebar-container").remove(),'
        'r._paper=r._paperdiv.insert("svg",":first-child").classed("main-svg",!0),r._toppaper=r._paperdiv.append("svg").classed("main-svg",!0),'
        'r._modebardiv=r._paperdiv.append("div"),delete r._modeBar,r._hoverpaper=r._paperdiv.append("svg").classed("main-svg",!0)')

    # check if text to be replaced is in current file
    if old not in filedata:
        # probable reason is other version of plotly/htmlwidgets
        logging.warning("Text in interactive plot could not be made selectable.\
                        Check versions of htmlwidgets and plotly")

    new = ('r._paperdiv=r._container.selectAll(".svg-container").data([0]),'
        'r._paperdiv.enter().append("div").classed("svg-container",!0).style("position","relative"),'
        'r._glcontainer=r._paperdiv.selectAll(".gl-container").data([{}]),r._glcontainer.enter().append("div").classed("gl-container",!0),'
        'r._paperdiv.selectAll(".main-svg").remove(),r._paperdiv.select(".modebar-container").remove(),'
        'r._paper=r._paperdiv.insert("svg",":first-child").classed("user-select-none",!0).classed("main-svg",!0),'
        'r._toppaper=r._paperdiv.append("svg").classed("user-select-none",!0).classed("main-svg",!0),'
        'r._modebardiv=r._paperdiv.append("div").classed("user-select-none",!0),delete r._modeBar,'
        'r._hoverpaper=r._paperdiv.append("svg").classed("user-select-none",!0).classed("main-svg",!0)')


    filedata = filedata.replace(old, new)

    with open(script_path, 'w') as file_out:
        file_out.write(filedata)


def make_selfcontained(html_file_path, html_dependencies_dir, html_out_path):
    """Make HTML file self-contained.

    Replaces links to scripts in the head part of HTML with
    script content to make it selfcontained

    Args:
      html_file_path(str): path to HTML plot file
      html_dependencies_dir(str): path to dir with HTML dependencies
      html_out_path(str): path to save self-contained HTML file to
    """

    with open(html_file_path) as html_file:
        soup = bs(html_file, 'html.parser')
        # find script elements in HTML file
        scripts = soup.find_all("script")
        for script_item in scripts:
            # check if element has src (link) attribute
            try:
                script_path = script_item.attrs['src']
            except:
                continue

            # script that is cruical to make text in hover selectable
            # and needs to be modified for that
            if 'plotly-main-' in script_path:
                make_selectable(html_dependencies_dir + script_path)

            # other scripts: read content and minify if it is javascript
            with open(html_dependencies_dir + script_path) as dependency_file:
                if script_path.endswith('.js'):
                    dependency_text = jsmin(dependency_file.read())
                else:
                    dependency_text = dependency_file.read()

            # remove script element from HTML
            script_item.extract()
            # add element with contend of script
            selfcontained_script = soup.new_tag('script')
            selfcontained_script.string = dependency_text
            soup.head.append(selfcontained_script)

        # replace links to css files with their content
        # elements are found by tag link
        head_links = soup.head.find_all('link')
        for link in head_links:
            # get path to css file
            css_path = link.attrs['href']

            # remove element from HTML
            link.extract()
            # retrieve content of css file
            with open(html_dependencies_dir + css_path) as css_file:
                styling_tag = soup.new_tag('style')
                styling_tag.string = css_file.read()

                # add content of css directly into HTML head
                soup.head.append(styling_tag)

    # write modified HTML to output location
    with open(html_out_path, 'w') as html_file:
        html_file.write(str(soup))


def change_title(output_path, html_out_path):
    """Change title of html file to species name.

    Reads species name from file temporarily created from
    taxonomic_assignment.py and replaces old title with species name

    Args:
      output_path: path to directory of MILTS report
      html_out_path: path to HTML plot file
    """

    # read species name from temporary file
    with open(output_path+'tmp/tmp.query_label', 'r') as tmp_file:
        next(tmp_file)
        query_name = next(tmp_file).strip()

    with open(html_out_path) as html_file:
        soup = bs(html_file, 'html.parser')
        # find and delete old title
        title = soup.find("title")
        title.extract()
        # replace with new title
        new_title = soup.new_tag('title')
        new_title.string = query_name
        soup.head.append(new_title)

    # write modified HTML to same location
    with open(html_out_path, 'w') as html_file:
        html_file.write(str(soup))


def perform_adjustments(cfg):
    """Collectively perform all modifications to HTML plot file.

    HTML is made self-contained, text is made selectable and
    title is changed to species name. Location of HTML is changed from
    tmp directory to final location in taxonomic_assignment dir

    Args:
      cfg: Config object with config parameters
    """

    # temporary location of HTML
    html_file_path = cfg.output_path + 'tmp/3D_plot.html'
    html_dependencies_dir = cfg.output_path + 'tmp/'
    # final location
    html_out_path = cfg.output_path + 'taxonomic_assignment/3D_plot.html'

    # makes selfcontained and calls script to make text selectable
    make_selfcontained(html_file_path, html_dependencies_dir, html_out_path)
    change_title(cfg.output_path, html_out_path)


def main():
    """Call module directly with preprocessed config file"""
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = prepare_and_check.cfg2obj(config_path)

    perform_adjustments(cfg)


if __name__ == '__main__':
    main()
