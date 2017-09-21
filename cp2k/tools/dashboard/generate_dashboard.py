#!/usr/bin/python
# -*- coding: utf-8 -*-

# Generates the CP2K Dashboard html page
# Inspired by Iain's cp2k_page_update.sh
#
# author: Ole Schuett

import sys
import os
import smtplib
from email.mime.text import MIMEText
import re
import gzip
from datetime import datetime, timedelta
import subprocess
import traceback
from os import path
from pprint import pformat
from xml.dom import minidom
from glob import glob
import itertools

try:
    from urllib.parse import urlencode
except ImportError:
    from urllib import urlencode

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

import matplotlib as mpl
mpl.use('Agg')  # change backend, to run without X11
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#===============================================================================
def main():
    if(len(sys.argv) not in (5, 6)):
        print("Usage update_dashboard.py <config-file> <addressbook> <status-file> <output-dir> [--full-archive]")
        sys.exit(1)

    config_fn, abook_fn, status_fn, outdir = sys.argv[1:5]
    assert(outdir.endswith("/"))
    assert(path.exists(config_fn))

    full_archive = False
    if(len(sys.argv) == 6):
        assert(sys.argv[5] == "--full-archive")
        full_archive = True

    config = configparser.ConfigParser()
    config.read(config_fn)

    if(full_archive):
        log = svn_log() # fetch entire history
        gen_archive(config, log, outdir, full_archive=True)
    else:
        log = svn_log(limit=100)
        gen_frontpage(config, log, abook_fn, status_fn, outdir)
        gen_archive(config, log, outdir)

#===============================================================================
def gen_frontpage(config, log, abook_fn, status_fn, outdir):
    addressbook = dict([line.split() for line in open(abook_fn).readlines()])

    if(path.exists(status_fn)):
        status = eval(open(status_fn).read())
    else:
        status = dict()

    trunk_revision = log[0]['num']
    log_index = dict([(r['num'], r) for r in log])
    now = datetime.utcnow().replace(microsecond=0)

    output  = html_header(title="CP2K Dashboard", rev=trunk_revision)
    output += '<div id="flex-container"><div>\n'
    output += html_svnbox(log)
    output += html_linkbox()
    output += '</div>\n'
    output += '<table border="1" cellspacing="3" cellpadding="5">\n'
    output += '<tr><th>Name</th><th>Host</th><th>Status</th>'
    output += '<th>Revision</th><th>Summary</th><th>Last OK</th><th>Tickets</th></tr>\n\n'

    def get_sortkey(s):
        return config.getint(s, "sortkey")

    for s in sorted(config.sections(), key=get_sortkey):
        print("Working on summary entry of: "+s)
        name        = config.get(s,"name")
        host        = config.get(s,"host")
        report_type = config.get(s,"report_type")
        report_url  = config.get(s,"report_url")
        do_notify   = config.getboolean(s,"notify") if(config.has_option(s,"notify")) else False
        timeout     = config.getint(s,"timeout") if(config.has_option(s,"timeout")) else 24

        # find latest revision that should have been tested by now
        freshness_threshold = now - timedelta(hours=timeout)
        revs_beyond_threshold = [ r['num'] for r in log if r['date'] < freshness_threshold ]
        threshold_rev = revs_beyond_threshold[0]

        # get and parse report
        report_txt = retrieve_report(report_url)
        report = parse_report(report_txt, report_type)

        if(s not in status.keys()):
            status[s] = {'last_ok': None, 'notified': False}

        if(report['status'] == "OK"):
            status[s]['last_ok'] = report['revision']
            status[s]['notified'] = False
        elif(do_notify and not status[s]['notified']):
            send_notification(report, addressbook, status[s]['last_ok'], log_index, name, s)
            status[s]['notified'] = True

        uptodate = report['revision']==trunk_revision # report from latest commit?
        if(report['revision'] != None):
            if(report['revision']<threshold_rev):
                report['status'] = "OUTDATED"
            elif(report['status'] in ("OK", "FAILED")):
                # store only useful and fresh reports, prevents overwriting archive
                fn = outdir+"archive/%s/rev_%d.txt.gz"%(s,report['revision'])
                write_file(fn, report_txt, gz=True)

        output += '<tr align="center">'
        output += '<td align="left"><a href="archive/%s/index.html">%s</a></td>'%(s, name)
        output += '<td align="left">%s</td>'%host
        output += status_cell(report['status'], report_url, uptodate)

        #Revision
        output += revision_cell(report['revision'], trunk_revision)

        #Summary
        output += '<td align="left">%s</td>'%report['summary']

        #Last OK
        if(report['status'] != "OK"):
            output += revision_cell(status[s]['last_ok'], trunk_revision)
        else:
            output += '<td></td>'

        output += ticket_cell(label=s)

        output += '</tr>\n\n'

    output += '</table>\n'
    output += '<div id="dummybox"></div></div>\n' # complete flex-container
    output += html_footer()
    write_file(outdir+"index.html", output.encode("utf8"))
    write_file(status_fn, pformat(status))

#===============================================================================
def gen_archive(config, log, outdir, full_archive=False):
    log_index = dict([(r['num'], r) for r in log])

    if(full_archive):
        print("Doing the full archive index pages")
        trunk_revision = None # trunk_version changes too quickly, leave it out.
        out_fn = "index_full.html"
        other_index_link = '<p>View <a href="index.html">recent archive</a></p>'
    else:
        print("Doing recent archive index pages")
        trunk_revision = log[0]['num']
        out_fn = "index.html"
        other_index_link = '<p>View <a href="index_full.html">full archive</a></p>'

    url_list = ""
    for s in config.sections():
        print("Working on archive page of: "+s)
        name        = config.get(s,"name")
        report_type = config.get(s,"report_type")
        info_url    = config.get(s,"info_url") if(config.has_option(s,"info_url")) else None

        # read all archived reports
        archive_reports = dict()
        for fn in sorted(glob(outdir+"archive/%s/rev_*.txt.gz"%s), reverse=True):
            report_txt = gzip.open(fn, 'rb').read()
            report = parse_report(report_txt, report_type)
            report['url'] = path.basename(fn)[:-3]
            archive_reports[report['revision']] = report

        # generate archive index
        archive_output = html_header(title=name, rev=trunk_revision)
        archive_output += '<p>Go back to <a href="../../index.html">main page</a></p>\n'
        if(info_url):
            archive_output += '<p>Get <a href="%s">more information</a></p>\n'%info_url
        if(report_type == "generic"):
            archive_output += gen_plots(archive_reports, log, outdir+"archive/"+s+"/", full_archive)
        archive_output += other_index_link
        archive_output += '<table border="1" cellspacing="3" cellpadding="5">\n'
        archive_output += '<tr><th>Revision</th><th>Status</th><th>Summary</th><th>Author</th><th>Commit Message</th></tr>\n\n'

        # loop over all relevant revisions
        rev_start = max(min(archive_reports.keys()), min(log_index.keys()))
        rev_end = max(log_index.keys())
        for r in range(rev_end, rev_start-1, -1):
            archive_output += '<tr>'
            archive_output += revision_cell(r, trunk_revision)
            if(archive_reports.has_key(r)):
                report = archive_reports[r]
                archive_output += status_cell(report['status'], report['url'])
                archive_output += '<td align="left">%s</td>'%report['summary']
                url_list += "https://dashboard.cp2k.org/archive/%s/%s.gz\n"%(s, report['url'])
            else:
                archive_output += 2*'<td></td>'
            svn_rev = log_index[r]
            archive_output += '<td align="left">%s</td>'%svn_rev['author']
            archive_output += '<td align="left">%s</td>'%svn_rev['msg'].split("\n")[0]
            archive_output += '</tr>\n\n'

        archive_output += '</table>\n'
        archive_output += other_index_link
        archive_output += html_footer()
        write_file(outdir+"archive/%s/%s"%(s,out_fn), archive_output.encode("utf8"))

    out_fn = "list_full.txt" if (full_archive) else "list_recent.txt"
    write_file(outdir+"archive/"+out_fn, url_list)

#===============================================================================
def gen_plots(all_reports, log, outdir, full_archive):
    # collect plot data
    plots = {}
    for revision in sorted(all_reports.keys()):
        report = all_reports[revision]
        for p in report['plots']:
            if(p['name'] not in plots.keys()):
                plots[p['name']] = {'curves':{}}
            plots[p['name']]['title'] = p['title'] # update title
            plots[p['name']]['ylabel'] = p['ylabel'] # update label
        for pp in report['plotpoints']:
            p = plots[pp['plot']]
            if(pp['name'] not in p['curves'].keys()):
                p['curves'][pp['name']] = {'x':[], 'y':[], 'yerr':[]}
            c = p['curves'][pp['name']]
            c['x'].append(report['revision'])
            c['y'].append(pp['y'])
            c['yerr'].append(pp['yerr'])
            c['label'] = pp['label'] # update label

    # write raw data
    tags = sorted([(pname, cname) for pname, p in plots.items() for cname in p['curves'].keys()])
    if(tags):
        raw_output = "#%9s"%"revision"
        for pname, cname in tags:
            raw_output += "   %18s   %22s"%(pname+"/"+cname,pname+"/"+cname+"_err")
        raw_output += "\n"
        for revision in sorted(all_reports.keys(), reverse=True):
            report = all_reports[revision]
            raw_output += "%10d"%revision
            for pname, cname in tags:
                pp = [pp for pp in report['plotpoints'] if(pp['plot']==pname and pp['name']==cname)]
                assert(len(pp)<=1)
                if(pp):
                    raw_output += "   %18f   %22f"%(pp[0]['y'],pp[0]['yerr'])
                else:
                    raw_output += "   %18s   %22s"%("?","?")
            raw_output += "\n"
        write_file(outdir+"plot_data.txt", raw_output)

    # create png images
    fig_ext = "_full.png" if(full_archive) else ".png"
    rev_end = log[0]['num']
    rev_start = min(all_reports.keys()) if(full_archive) else rev_end-100
    for pname, p in plots.items():
        print("Working on plot: "+pname)
        fig = plt.figure(figsize=(12,4))
        fig.subplots_adjust(bottom=0.18, left=0.06, right=0.70)
        fig.suptitle(p['title'], fontsize=14, fontweight='bold', x=0.4)
        ax = fig.add_subplot(111)
        ax.set_xlabel('SVN Revision')
        ax.set_ylabel(p['ylabel'])
        markers = itertools.cycle('os>^*')
        for cname in sorted(p['curves'].keys()):
            c= p['curves'][cname]
            if(full_archive):
                ax.plot(c['x'], c['y'], label=c['label'], linewidth=2) # less crowded
            else:
                ax.errorbar(c['x'], c['y'], yerr=c['yerr'], label=c['label'],
                            marker=markers.next(), linewidth=2, markersize=6)
        ax.set_xlim(rev_start-1, rev_end+1)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left',
                  numpoints=1, fancybox=True, shadow=True, borderaxespad=0.0)
        visibles = [[y for x,y in zip(c['x'],c['y']) if x>=rev_start] for c in p['curves'].values()] # visible y-values
        ymin  = min([min(ys) for ys in visibles]) # lowest point from lowest curve
        ymax = max([max(ys) for ys in visibles]) # highest point from highest curve
        if(full_archive):
            ax.set_ylim(0.98*ymin, 1.02*ymax)
        else:
            ymax2 = max([min(ys) for ys in visibles]) # lowest point from highest curve
            ax.set_ylim(0.98*ymin, min(1.02*ymax, 1.3*ymax2))  # protect against outlayers
        fig.savefig(outdir+pname+fig_ext)

    # write html output
    html_output = ""
    for pname in sorted(plots.keys()):
        html_output += '<a href="plot_data.txt"><img src="%s" alt="%s"></a>\n'%(pname+fig_ext, plots[pname]['title'])
    return(html_output)

#===============================================================================
def retrieve_report(report_url):
    try:
        return urlopen(report_url, timeout=5).read()
    except:
        print(traceback.print_exc())
        return None

#===============================================================================
def parse_report(report_txt, report_type):
    if(report_txt==None):
        return( {'status':'UNKNOWN', 'summary':'Error while retrieving report.', 'revision':None} )
    try:
        if(report_type == "regtest"):
            return parse_regtest_report(report_txt)
        elif(report_type == "generic"):
            return parse_generic_report(report_txt)
        else:
            raise(Exception("Unknown report_type"))
    except:
        print(traceback.print_exc())
        return( {'status':'UNKNOWN', 'summary':'Error while parsing report.', 'revision':None} )

#===============================================================================
def send_notification(report, addressbook, last_ok, svn_log, name, s):
    rev_end = report['revision'] if(report['revision']) else max(svn_log.keys())
    if(rev_end==last_ok): return # probably a flapping tester
    authors = set([svn_log[r]['author'] for r in range(last_ok+1, rev_end+1)])
    emails = [addressbook[a] for a in authors]
    print("Sending email to: "+", ".join(emails))

    msg_txt  = "Dear CP2K developer,\n\n"
    msg_txt += "the dashboard has detected a problem that one of your recent commits might have introduced.\n\n"
    msg_txt += "   test name:      %s\n"%name
    msg_txt += "   report state:   %s\n"%report['status']
    msg_txt += "   report summary: %s\n"%report['summary']
    msg_txt += "   last OK rev:    %d\n\n"%last_ok
    msg_txt += "For more information visit:\n"
    msg_txt += "   https://dashboard.cp2k.org/archive/%s/index.html \n\n"%s
    msg_txt += "Sincerely,\n"
    msg_txt += "  your CP2K Dashboard ;-)\n"

    msg = MIMEText(msg_txt)
    msg['Subject'] = "Problem with "+name
    msg['From']    = "CP2K Dashboard <dashboard@cp2k.org>"
    msg['To']      = ", ".join(emails)

    smtp_conn = smtplib.SMTP('localhost')
    smtp_conn.sendmail(msg['From'], emails, msg.as_string())
    smtp_conn.quit()

#===============================================================================
def html_header(title, rev=None):
    output  = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n'
    output += '<html><head>\n'
    output += '<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n'
    output += '<meta http-equiv="refresh" content="200">\n'
    output += '<link rel="icon" type="image/x-icon" href="data:image/x-icon;base64,AAABAAEAEBAQAAAAAAAoAQAAFgAAACgAAAAQAAAAIAAAAAEABAAAAAAAgAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAD/AAmRCQAAb/8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAiIgAAAAAAACIiAAAAAAAAIiIAAAAAAAAAAAAAAAAAADMzAAAAAAAAMzMAAAAAAAAzMwAAAAAAAAAAAAAAAAAAEREAAAAAAAAREQAAAAAAABERAAAAAAAAAAAAAAD+fwAA/n8AAPw/AAD4HwAA+B8AAPgfAAD4HwAA+B8AAPgfAAD4HwAA+B8AAPgfAAD4HwAA+B8AAPgfAAD4HwAA">\n'
    output += '<style type="text/css">\n'
    output += '.ribbon {\n'
    output += '  overflow: hidden;\n'
    output += '  position: absolute;\n'
    output += '  right:0px;\n'
    output += '  top: 0px;\n'
    output += '  width: 200px;\n'
    output += '  height: 200px;\n'
    output += '}\n'
    output += '.ribbon a {\n'
    output += '  position: relative;\n'
    output += '  white-space: nowrap;\n'
    output += '  background-color: #a00;\n'
    output += '  border: 1px solid #faa;\n'
    output += '  color: #fff;\n'
    output += '  display: block;\n'
    output += '  font: bold 11pt sans-serif;\n'
    output += '  padding: 7px;\n'
    output += '  top: 35px;\n'
    output += '  right: 10px;\n'
    output += '  width: 300px;\n'
    output += '  text-align: center;\n'
    output += '  text-decoration: none;\n'
    output += '  transform: rotate(45deg);\n'
    output += '  box-shadow: 0 0 10px #888;\n'
    output += '}\n'
    output += '#flex-container {\n'
    output += '  display: -webkit-flex; /* Safari */\n'
    output += '  display: flex;\n'
    output += '  -webkit-flex-flow: row wrap-reverse; /* Safari */\n'
    output += '  flex-flow:         row wrap-reverse;\n'
    output += '  -webkit-justify-content: space-around; /* Safari */\n'
    output += '  justify-content:         space-around;\n'
    output += '  -webkit-align-items: flex-end; /* Safari */\n'
    output += '  align-items:         flex-end;\n'
    output += '}\n'
    output += '.sidebox {\n'
    output += '  width: 15em;\n'
    output += '  border-radius: 1em;\n'
    output += '  box-shadow: .2em .2em .7em 0 #777;\n'
    output += '  background: #f7f7f0;\n'
    output += '  padding: 1em;\n'
    output += '  margin: 40px 20px;\n'
    output += '}\n'
    output += '.sidebox h2 {\n'
    output += '  margin: 0 0 0.5em 0;\n'
    output += '}\n'
    output += '.sidebox p {\n'
    output += '  margin: 0.5em;\n'
    output += '}\n'
    output += '#dummybox {\n'
    output += '  width: 15em;\n'
    output += '}\n'
    output += '</style>\n'
    output += '<title>%s%s</title>\n'%(title, (" (%d)"%rev if rev else ""))
    output += '</head><body>\n'
    output += '<div class="ribbon"><a href="https://cp2k.org/dev:dashboard">Need Help?</a></div>\n'
    output += '<center><h1>%s</h1></center>\n'%title.upper()
    return(output)

#===============================================================================
def html_linkbox():
    output  = '<div class="sidebox">\n'
    output += '<h2>More...</h2>\n'
    output += '<a href="regtest_survey.html">Regtest Survey</a><br>\n'
    output += '<a href="https://www.cp2k.org/static/coverage/">Test Coverage</a><br>\n'
    output += '<a href="discontinued_tests.html">Discontinued Tests</a><br>\n'
    output += '</div>\n'
    return(output)

#===============================================================================
def html_svnbox(log):
    now = datetime.utcnow()
    output  = '<div class="sidebox">\n'
    output += '<h2>Recent Commits</h2>\n'
    for r in log[0:10]:
        url = "https://sourceforge.net/p/cp2k/code/%d/"%r['num']
        msg = r['msg'].split("\n")[0]
        if(len(msg) > 27):
            msg = msg[:26] + "..."
        output += '<p><a title="%s" href="%s">%s</a><br>\n'%(r['msg'], url, msg)
        delta = now - r['date']
        age = delta.days*24.0 + delta.seconds/3600.0
        output += '<small>%s %s %.1fh ago.</small></p>\n'%(r['num'], r['author'], age)
    output += '</div>\n'
    return(output)

#===============================================================================
def html_footer():
    now = datetime.utcnow().replace(microsecond=0)
    output  = '<p><small>Page last updated: %s</small></p>\n'%now.isoformat()
    output += '</body></html>'
    return(output)

#===============================================================================
def write_file(fn, content, gz=False):
    d = path.dirname(fn)
    if(len(d) > 0 and not path.exists(d)):
        os.makedirs(d)
        print("Created dir: "+d)
    f = gzip.open(fn, 'wb') if(gz) else open(fn, "w")
    f.write(content)
    f.close()
    print("Wrote: "+fn)

#===============================================================================
def svn_log(limit=None):
    sys.stdout.write("Fetching svn log... ")
    sys.stdout.flush()
    # xml version contains nice UTC timestamp
    if(limit):
        cmd = "svn log --limit %d svn://svn.code.sf.net/p/cp2k/code --xml"%limit
    else:
        # our server is much faster, but might be 5 minutes old
        cmd = "svn log https://svn.cp2k.org/cp2k --xml"
    log_xml = check_output(cmd.split())
    dom = minidom.parseString(log_xml)
    revisions = []
    for entry in dom.getElementsByTagName("logentry"):
        rev = dict()
        rev['num'] = int(entry.attributes['revision'].value)
        rev_date_str = entry.getElementsByTagName('date')[0].firstChild.nodeValue
        rev['date'] = datetime.strptime(rev_date_str[:19], '%Y-%m-%dT%H:%M:%S')
        # some revisions lack an author
        author = entry.getElementsByTagName('author')
        rev['author'] = author[0].firstChild.nodeValue if (len(author)==1) else ""
        msg = entry.getElementsByTagName('msg')[0].firstChild
        rev['msg'] = msg.nodeValue.strip() if(msg) else ""
        revisions.append(rev)
    print("done.")
    return(revisions)

#===============================================================================
def status_cell(status, report_url, uptodate=True):
    if(status == "OK"):
        bgcolor = "#00FF00" if(uptodate) else "#8CE18C"
    elif(status == "FAILED"):
        bgcolor = "#FF0000" if(uptodate) else "#E18C8C"
    else:
        bgcolor = "#d3d3d3"
    return('<td bgcolor="%s"><a href="%s">%s</a></td>'%(bgcolor, report_url, status))

#===============================================================================
def revision_cell(rev, trunk_rev):
    if(rev == None):
        return('<td>N/A</td>')
    rev_url = "https://sourceforge.net/p/cp2k/code/%d/"%rev
    rev_delta = "(%d)"%(rev - trunk_rev) if(trunk_rev) else ""
    output = '<td align="left"><a href="%s">%s</a> %s</td>'%(rev_url, rev, rev_delta)
    return(output)

#===============================================================================
def ticket_cell(label):
    base_url = "https://sourceforge.net/p/cp2k/bugs"
    new_url = base_url+"/new/?" + urlencode({'labels':label})
    query = urlencode({'q':'!status:wont-fix && !status:closed && labels:"%s"'%label})
    feed_url = base_url+"/search_feed/?limit=25&sort=ticket_num_i+asc&" + query
    output = '<td  align="right">'
    try:
        # sometime the http-request to sourceforge times out
        tickets_xml = urlopen(feed_url, timeout=5).read()
        dom = minidom.parseString(tickets_xml)
        for entry in dom.getElementsByTagName("item"):
            title = entry.getElementsByTagName('title')[0].firstChild.nodeValue
            link = entry.getElementsByTagName('link')[0].firstChild.nodeValue
            tid = int(link.strip("/ ").split("/")[-1])
            output += '<a href="%s" title="%s">#%d</a>, '%(link, title, tid)
    except:
        print(traceback.print_exc())
        output += "N/A "
    output += '<a href="%s"'%new_url
    output += ' style="text-decoration:none;font-weight:bold;font-size:larger;"'
    output += ' title="Create a new Ticket">+</a></td>'
    return(output)

#===============================================================================
def parse_regtest_report(report_txt):
    m = re.search("svn: .*(Can't connect .*):", report_txt)
    if(m):
        return({'revision':None, 'status':'UNKNOWN', 'summary':m.group(1)})

    report = dict()
    report['revision'] = int(re.search("(revision|Revision:) (\d+)\.?\n", report_txt).group(2))

    if("LOCKFILE" in report_txt):
        report['status'] = "UNKNOWN"
        report['summary'] = "Test directory is locked."
        return(report)

    m = re.search("\nGREPME (\d+) (\d+) (\d+) (\d+) (\d+) (.+)\n", report_txt)
    if(not m and re.search("make: .* Error .*", report_txt)):
        report['status'] = "FAILED"
        report['summary'] = "Compilation failed."
        return(report)

    runtime_errors = int(m.group(1))
    wrong_results  = int(m.group(2))
    correct_tests  = int(m.group(3))
    new_inputs     = int(m.group(4))
    num_tests      = int(m.group(5))
    memory_leaks   = int(m.group(6).replace("X", "0"))

    report['summary'] = "correct: %d / %d"%(correct_tests, num_tests)
    if(new_inputs > 0):
        report['summary'] += "; new: %d"%new_inputs
    if(wrong_results > 0):
        report['summary'] += "; wrong: %d"%wrong_results
    if(runtime_errors > 0):
        report['summary'] += "; failed: %d"%runtime_errors
    if(memory_leaks > 0):
        report['summary'] += "; memleaks: %d"%memory_leaks

    runtimes = [float(m) for m in re.findall("\nRegtest took (.+) seconds.\n", report_txt)]
    report['summary'] += "; %.0fmin"%(sum(runtimes)/60.0)

    if(wrong_results>0 or runtime_errors>0 or memory_leaks>0):
        report['status'] = "FAILED"
    else:
        report['status'] = "OK"

    return(report)

#===============================================================================
def parse_generic_report(report_txt):
    m = re.search("svn: .*(Can't connect .*):", report_txt)
    if(m):
        return({'revision':None, 'status':'UNKNOWN', 'summary':m.group(1)})
    report = dict()
    report['revision']  = int(re.search("(^|\n)Revision: (\d+)\n", report_txt).group(2))
    report['summary'] = re.search("(^|\n)Summary: (.+)\n", report_txt).group(2)
    report['status'] = re.search("(^|\n)Status: (.+)\n", report_txt).group(2)
    report['plots'] = [eval("dict(%s)"%m[1]) for m in re.findall("(^|\n)Plot: (.+)(?=\n)", report_txt)]
    report['plotpoints'] = [eval("dict(%s)"%m[1]) for m in re.findall("(^|\n)PlotPoint: (.+)(?=\n)", report_txt)]
    return(report)

#===============================================================================
def check_output(command, **kwargs):
    p = subprocess.Popen(command, stdout=subprocess.PIPE, **kwargs)
    output = p.communicate()[0]
    assert(p.wait() == 0)
    return(output)

#===============================================================================
if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
    pass #TODO implement selftest
else:
    main()
#EOF
