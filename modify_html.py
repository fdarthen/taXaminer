# -*- coding: utf-8 -*-

from bs4 import BeautifulSoup as bs
from jsmin import jsmin
import sys
import yaml



def make_selectable(script_path):
    """ this is a highly specified function and may not work with other versions of htmlwidgets and plotly
     (developed for htmlwidgets-1.5.7 and plotly-4.9.0; compatiblity tested with htmlwidgets-1.5.4 and plotly-4.10.0)

     replaces parts of javascript script to make the text in the hoverwindow selectable """

    with open(script_path, 'r') as file_in:
        filedata = file_in.read()

    old = ('r._paperdiv=r._container.selectAll(".svg-container").data([0]),'
        'r._paperdiv.enter().append("div").classed("user-select-none",!0).classed("svg-container",!0).style("position","relative"),'
        'r._glcontainer=r._paperdiv.selectAll(".gl-container").data([{}]),r._glcontainer.enter().append("div").classed("gl-container",!0),'
        'r._paperdiv.selectAll(".main-svg").remove(),r._paperdiv.select(".modebar-container").remove(),'
        'r._paper=r._paperdiv.insert("svg",":first-child").classed("main-svg",!0),r._toppaper=r._paperdiv.append("svg").classed("main-svg",!0),'
        'r._modebardiv=r._paperdiv.append("div"),delete r._modeBar,r._hoverpaper=r._paperdiv.append("svg").classed("main-svg",!0)')

    if old not in filedata:
        print("Error: text in interactive plot could not be made selectable.")
        # probable reason is new version of plotly/htmlwidgets

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


def make_selfcontained(html_file_path, html_dependencies_grounddir, html_out_path):
    """ replace links to scripts in head of html with script content to make selfcontained """

    with open(html_file_path) as html_file:
        soup = bs(html_file, 'html.parser')
        scripts = soup.find_all("script")
        for script_item in scripts:
            try:
                script_path = script_item.attrs['src']
            except:
                continue

            if 'plotly-main-' in script_path:
                make_selectable(html_dependencies_grounddir + script_path)


            with open(html_dependencies_grounddir + script_path) as dependency_file:
                if script_path.endswith('.js'):
                    dependency_text = jsmin(dependency_file.read())
                else:
                    dependency_text = dependency_file.read()

            script_item.extract()
            selfcontained_script = soup.new_tag('script')
            selfcontained_script.string = dependency_text
            soup.head.append(selfcontained_script)


        head_links = soup.head.find_all('link')
        for link in head_links:
            css_path = link.attrs['href']

            link.extract()
            with open(html_dependencies_grounddir + css_path) as css_file:
                styling_tag = soup.new_tag('style')
                styling_tag.string = css_file.read()

                soup.head.append(styling_tag)

    with open(html_out_path, 'w') as html_file:
        html_file.write(str(soup))


def change_title(output_path, html_out_path):
    """ change title of html file to species name """

    with open(output_path+'tmp/tmp.query_label', 'r') as tmp_file:
        next(tmp_file)
        query_name = next(tmp_file).strip()

    with open(html_out_path) as html_file:
        soup = bs(html_file, 'html.parser')
        title = soup.find("title")
        title.extract()
        new_title = soup.new_tag('title')
        new_title.string = query_name
        soup.head.append(new_title)

    with open(html_out_path, 'w') as html_file:
        html_file.write(str(soup))


def main():

    config_path = sys.argv[1]

    # read parameters from config file
    config_obj=yaml.safe_load(open(config_path,'r'))
    output_path=config_obj['output_path'] # complete output path (ENDING ON A SLASH!)

    html_file_path = output_path + 'tmp/3D_plot.html'
    html_dependencies_grounddir = output_path + 'tmp/'
    html_out_path = output_path + 'taxonomic_assignment/3D_plot.html'

    make_selfcontained(html_file_path, html_dependencies_grounddir, html_out_path)
    change_title(output_path, html_out_path)


if __name__ == '__main__':
    main()
