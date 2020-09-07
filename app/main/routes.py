from datetime import datetime
from flask import render_template, flash, redirect, url_for, request, g, jsonify, current_app
from flask_login import current_user, login_required
from flask_babel import _, get_locale
from flask.helpers import get_root_path
from guess_language import guess_language
import app
from app import db
from app.main.forms import EditProfileForm, EmptyForm, PostForm , SearchForm, \
                           VirusForm, MutationMapForm, LabResultsForm
from app.models import User, Post
from app.translate import translate
from app.main import bp
import dash

from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from operator import itemgetter #, attrgetter
import copy, os, random, sys, re, math
from os.path import dirname, join
import numpy as np
import pandas as pd
from collections import OrderedDict

from bokeh.io      import curdoc # , output_notebook, show
from bokeh.layouts import column, layout

from bokeh.core.properties import value
from bokeh.io       import output_file
from bokeh.plotting import figure
from bokeh.models   import ColumnDataSource, LabelSet, HoverTool, Grid, LinearAxis, Plot, VBar, Label
from bokeh.models   import Div, Select, Slider, TextInput, CustomJS, MultiSelect, Button, AutocompleteInput
from bokeh.models   import Range1d, FixedTicker, CheckboxGroup
from bokeh.events   import ButtonClick

sys.path.insert(1, '../src/')
from Basic          import *
from Sequence       import *
from Shannon        import *
from covid_seqs_lib import *

# source must return 2 dirs
root0        = "../../colaboracoes/covid/"
root_data0   = root0 + "fasta/sarscov2_202007/"
root_data    = os.path.join(root_data0, 'msa_0713_202007/')
root_protein = os.path.join(root_data, "protein/")
root_country = os.path.join(root_data, "countries/")

root_cut     = os.path.join(root_data, "cut/")
root_entropy = os.path.join(root_cut, "entropy/")
root_figure  = os.path.join(root_cut, 'figure/')
root_html    = os.path.join(root_cut, "html/")

root_templates = os.path.join(root0, "templates/")
root_ensembl   = os.path.join(root0, 'fasta/ensembl')
root_ncbi      = os.path.join(root0, 'fasta/ncbi')

root_protent_figure  = os.path.join(root_protein, 'figure/')
root_protent_entropy = os.path.join(root_protein, 'entropy/')
root_protent_html    = os.path.join(root_protein, "html/")

try:
    os.mkdir(root_protent_figure)
except:
    pass

try:
    os.mkdir(root_protent_html)
except:
    pass

filefas  = 'msa_sarscov2_0713.fasta'
filehead = 'msa_countries_ids_edited.tsv'

prjName = "GISAID"
typeSeq = "SARS_COV2"
pattern = ""

isProtein = True; isAligned = True; isFiltered = True

filesummary = 'summary.tsv'
filemeta    = 'metadata.tsv'
#--- Plot entropy barplot with Botek
base_filename = 'msa_sarscov2_202007_%s_protein_%s.fasta'
base_pickle   = 'msa_sarscov2_202007_%s_protein_%s.pickle'
xlabel = "AA positions"; moltype='Protein'
ylabel = 'H (bits)'


@bp.before_app_request
def before_request():
    if current_user.is_authenticated:
        current_user.last_seen = datetime.utcnow()
        db.session.commit()
        g.search_form = SearchForm()
    g.locale = str(get_locale())


@bp.route('/', methods=['GET', 'POST'])
@bp.route('/index', methods=['GET', 'POST'])
@login_required
def index():
    form = PostForm()
    if form.validate_on_submit():
        language = guess_language(form.post.data)
        if language == 'UNKNOWN' or len(language) > 5:
            language = ''
        post = Post(body=form.post.data, author=current_user,
                    language=language)
        db.session.add(post)
        db.session.commit()
        flash(_('Your post is now live!'))
        return redirect(url_for('main.index'))
    page = request.args.get('page', 1, type=int)
    posts = current_user.followed_posts().paginate(
        page, current_app.config['POSTS_PER_PAGE'], False)
    next_url = url_for('main.index', page=posts.next_num) \
        if posts.has_next else None
    prev_url = url_for('main.index', page=posts.prev_num) \
        if posts.has_prev else None

    print("####", current_app.file_index_css)
    return render_template('index.html', title=_('Home'), form=form,
                           posts=posts.items, next_url=next_url,
                           prev_url=prev_url, file_index_css=current_app.file_index_css)


@bp.route('/find_css_file')
def find_css_file():
    _root = 'app/static/css/'
    print(">>>> <<< ", os.path.join(_root, current_app.file_index_css))
    return os.path.join(_root, current_app.file_index_css)

@bp.route('/explore')
@login_required
def explore():
    page = request.args.get('page', 1, type=int)
    posts = Post.query.order_by(Post.timestamp.desc()).paginate(
        page, current_app.config['POSTS_PER_PAGE'], False)
    next_url = url_for('main.explore', page=posts.next_num) \
        if posts.has_next else None
    prev_url = url_for('main.explore', page=posts.prev_num) \
        if posts.has_prev else None
    return render_template('index.html', title=_('Explore'),
                           posts=posts.items, next_url=next_url,
                           prev_url=prev_url)

@bp.route('/user/<username>')
@login_required
def user(username):
    user = User.query.filter_by(username=username).first_or_404()
    page = request.args.get('page', 1, type=int)
    posts = user.posts.order_by(Post.timestamp.desc()).paginate(
        page, current_app.config['POSTS_PER_PAGE'], False)
    next_url = url_for('main.user', username=user.username,
                       page=posts.next_num) if posts.has_next else None
    prev_url = url_for('main.user', username=user.username,
                       page=posts.prev_num) if posts.has_prev else None
    form = EmptyForm()
    return render_template('user.html', user=user, posts=posts.items,
                           next_url=next_url, prev_url=prev_url, form=form)


@bp.route('/edit_profile', methods=['GET', 'POST'])
@login_required
def edit_profile():
    form = EditProfileForm(current_user.username)
    if form.validate_on_submit():
        current_user.username = form.username.data
        current_user.about_me = form.about_me.data
        db.session.commit()
        flash(_('Your changes have been saved.'))
        return redirect(url_for('main.edit_profile'))
    elif request.method == 'GET':
        form.username.data = current_user.username
        form.about_me.data = current_user.about_me
    return render_template('edit_profile.html', title=_('Edit Profile'), form=form)

@bp.route('/lab_results_202008')
@login_required
def lab_results_202008():
    experiments = ['Chondrocytes (Paola)', "Synoviocites (Isadora)", "Osteoclast (Eduardo)"]
    return render_template('lab_results_202008.html', title='Screening results', experiments=experiments)

@bp.route('/show_exp_results/<string:experiment>')
@login_required
def show_exp_results(experiment):
    print(">>> experiment", experiment)

    # register_dashapps()
    print(url_for('main.index'))
    print(url_for('/dashboard/'))

    return redirect(url_for('/dashboard/') )

    # return render_template('heatmap.html', title='Results for '+experiment, experiment=experiment)


# https://medium.com/@olegkomarov_77860/how-to-embed-a-dash-app-into-an-existing-flask-app-ea05d7a2210b
def register_dashapps():
    from app.dashPlots.layout    import layout
    from app.dashPlots.callbacks import register_callbacks

    # Meta tags for viewport responsiveness
    meta_viewport = {"name": "viewport", "content": "width=device-width, initial-scale=1, shrink-to-fit=no"}

    mydash = dash.Dash(__name__,
                         server=app,
                         url_base_pathname='/dashboard/',
                         assets_folder=get_root_path(__name__) + '/dashboard/assets/',
                         meta_tags=[meta_viewport])

    with app.app_context():
        mydash.title = 'Dashapp 1'
        mydash.layout = layout
        register_callbacks(mydash)

    _protect_dashviews(mydash)

def _protect_dashviews(dashapp):
    for view_func in dashapp.server.view_functions:
        if view_func.startswith(dashapp.config.url_base_pathname):
            dashapp.server.view_functions[view_func] = login_required(dashapp.server.view_functions[view_func])


@bp.route('/virus')
@login_required
def virus():
    form = VirusForm()
    return render_template('virus_sarscov2.html', title=_('SARS-CoV-2'), form=form)

@bp.route('/referencias_links')
def referencias_links():
    return render_template('index.html', title=_('SARS-CoV-2'), form=form)

@bp.route('/mutation_map')
def mutation_map():
    form = MutationMapForm()

    # dfh, countries = split_gisaid_fasta_headers(filehead, filefas, root_data, force=False)
    dfsumm = pdreadcsv(filesummary, root_protein, verbose=False)
    countries = list(dfsumm.country.unique())

    #--- metadata file
    dfmeta = pdreadcsv(filemeta, root_data0, verbose=False)


    #------------ objects ----------------------------------
    deschtml = Div(text=open("app/templates/mutmap.html").read(), sizing_mode="stretch_width")

    dicdata0 = {'country' : [], 'x' : [], 'y': [], 'aas': [], 'nns': []}
    dicdata1 = {'country' : [], 'x' : [], 'y': [], 'aas': [], 'nns': []}
    dicdata2 = {'country' : [], 'x' : [], 'y': [], 'aas': [], 'nns': []}

    source0 = ColumnDataSource(data=dicdata0)
    source1 = ColumnDataSource(data=dicdata1)
    source2 = ColumnDataSource(data=dicdata2)
    sources = [source0, source1, source2]

    increment = -1; ylabel = 'H (bits)'
    delwidth = .25

    colors = ['blue', 'red', 'green']

    key = 'S'
    countries2 = ['a', 'b', 'c' ]
    title0 = "Entropy for SARS-CoV-2 protein '%s'"
    title = title0%(key)

    TOOLTIPS=[
        ("country", "@country"),
        ("x", '@x{int}' ),
        ("aa", '@aas' ),
        ("#", "@nns"),
        ("h", '@y')
    ]

    #   min_border=0, toolbar_location=None, sizing_mode="scale_both"
    plot_width=1000; plot_height=600
    p = figure(title=title, plot_width=plot_width, plot_height=plot_height, tooltips=TOOLTIPS, sizing_mode="scale_both")

    # handles = []
    for i in range(3):
        handle = p.vbar(x='x', top='y', width=.25, color=colors[i], source=sources[i] ) # , legend=countries2[i]
        p.add_tools(HoverTool(renderers=[handle], tooltips=TOOLTIPS,  mode='vline')  )
        # handles.append(handle)


    p.x_range.range_padding = 0.6
    p.xgrid.grid_line_color = None
    # p.legend.location = "top_left"
    # p.legend.orientation = "vertical"

    xaxis = LinearAxis()
    # p.add_layout(xaxis, 'below') # major_label_overrides = xdic

    yaxis = LinearAxis()
    # p.add_layout(yaxis, 'left')

    p.add_layout(Grid(dimension=0, ticker=xaxis.ticker))
    p.add_layout(Grid(dimension=1, ticker=yaxis.ticker))

    p.xaxis.axis_label = xlabel
    p.yaxis.axis_label = ylabel

    # p.xaxis.major_label_overrides = xdic
    p.xaxis.major_label_orientation = math.pi/2

    # p.x_range=Range1d(0, xmaxi+5)
    # p.y_range=Range1d(0, maximus*1.2)

    p.title.text_font_size = '20pt'
    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"

    #---------------- controls ---------------------------------------
    # source_countries = ColumnDataSource(data={'country': countries})
    # def showmessage(title, content):
    #    dialogs(title, content, closable=True)

    txt_date_inf = TextInput(value="01/04/2020", title="Inferior (inclusive):")
    txt_date_sup = TextInput(value="12/05/2020", title="Superior (exclusive):")
    txt_dummy    = TextInput(value="---", title="dummy")
    txt_alert    = TextInput(value="---", title="alert")

    # min_month = Slider(title="initial month (2020)", start=1, end=6, value=4, step=1)
    # max_month = Slider(title="end month (2020)", start=1, end=7, value=6, step=1)

    '''
    checkbox_countries = CheckboxGroup(labels=countries, active=[0, 1, 2])
    '''
    autocomp0 = AutocompleteInput(completions=countries, value='Brazil',    background=colors[0])
    autocomp1 = AutocompleteInput(completions=countries, value='Argentina', background=colors[1])
    autocomp2 = AutocompleteInput(completions=countries, value='Chile',     background=colors[2])


    label_country0 = Div(text="""<br><b><span style="color: blue; ">Chose country:</span></b>""")
    label_country1 = Div(text="""    <b><span style="color: red;  ">Chose country:</span></b>""")
    label_country2 = Div(text="""    <b><span style="color: green;">Chose country:</span></b>""")
    group_countries = column(label_country0, autocomp0, label_country1, autocomp1, label_country2, autocomp2)

    '''
    reset_all = Button(label="reset all")

    def clear_checkbox_countries():
        checkbox_countries.active = [] * len(countries)

    reset_all.on_click(clear_checkbox_countries)

    group_countries = column(reset_all, checkbox_countries)
    '''
    keysall = keys + nspkeys
    # n = [i for i in range(len(keysall)) if keysall[i] = 'S'][0]
    select_protein = Select(options=keysall, value='S')

    label_protein = Div(text="""<br><b><span style="color: blue;">Chose a Protein/NSP</span></b><br>""")
    group_protein = column(label_protein, select_protein)

    label_date  = Div(text="""<b><span style="color: blue;">Chose a date interval</span></b><br>""")
    group_dates = column(label_date, txt_date_inf, txt_date_sup)

    button_submit = Button(label="refresh") # label='refresh',button_type='success'

    #---------------- udapte functions -------------------------------
    def update():

        _date_inf, _date_sup = select_ranges()
        if _date_inf is None: return

        print(_date_inf, _date_sup)

        global countries2
        '''
        countries2 = [checkbox_countries.labels[i] for i in checkbox_countries.active]
        if len(countries2) > 3:
            countries2 = countries2[:3]'''
        countries2 = [autocomp0.value, autocomp1.value, autocomp2.value]

        print(",".join(countries2))

        global key
        key = select_protein.value

        print("calc0", countries2[0])
        dicdata, nSamples0, maxVal0 = select_samples_by_country(countries2[0], 0, _date_inf, _date_sup)
        source0.data = dicdata

        print("calc1", countries2[1])
        dicdata, nSamples1, maxVal1 = select_samples_by_country(countries2[1], 1, _date_inf, _date_sup)
        source1.data = dicdata

        print("calc2", countries2[2])
        dicdata, nSamples2, maxVal2 = select_samples_by_country(countries2[2], 2, _date_inf, _date_sup)
        source2.data = dicdata

        maxVal = np.max( [maxVal0, maxVal1, maxVal2] )

        # .strftime('%d/%b/%Y')
        mat = [x + "(%d)" for x in countries2]
        stri = ", ".join(mat)%(nSamples0, nSamples1, nSamples2)

        p.title.text = "Protein '%s', countries: %s between %s and %s"%(key, stri, _date_inf, _date_sup )
        curdoc().title = "Protein '%s', countries: %s"%(key, ", ".join(countries2))


        '''
        try:
            del(hcirc0); del(hcirc1); del(hcirc2)
        except:
            pass
        '''

        '''
        hcirc0 = p.circle( x=0, y=0, radius=0, color=colors[0], legend_label= countries2[0])
        hcirc1 = p.circle( x=0, y=0, radius=0, color=colors[1], legend_label= countries2[1])
        hcirc2 = p.circle( x=0, y=0, radius=0, color=colors[2], legend_label= countries2[2])

        p.legend.location    = "top_left"
        p.legend.orientation = "vertical"
        '''

    '''
    def update_graph():
        p.legend.items = [ (countries2[0], handles[0]),
                           (countries2[1], handles[1]),
                           (countries2[2], handles[2])  ]

    '''

    def select_ranges():
        _date_inf = txt_date_inf.value
        _date_sup = txt_date_sup.value
        if len(_date_inf) != 10 or len(_date_sup) != 10:
            title = 'Format date'
            content = 'Date must be in the following format: dd/mm/yyyy'
            return None, None

        '''
        _date_inf = '01/%d/2020'%(min_month.value)
        _date_sup = '01/%d/2020'%(max_month.value)
        '''
        return _date_inf, _date_sup


    def select_samples_by_country(country, i, _date_inf, _date_sup, ini = 0, isProtein=True, increment = -1):
        keyName = key

        shan = Shannon(prjName, isProtein, country, keyName, base_filename,
                       root_data, root_protein, root_protent_figure, root_protent_entropy, root_protent_html, root_templates)

        #--- load metadata object
        dic = pdloaddic(shan.file_entropy, path=root_protent_entropy, verbose=False)
        # print(dic.keys())
        df_metadata = dic['metadata']

        fileprot = base_filename%(country, key)
        filefull = os.path.join(root_protein, fileprot)

        #--- load protein object
        mseqProt = MySequence(prjName, root_protein)
        mseqProt.readFasta(filefull, showmessage=False)

        #-- filter de i's related to the filter - convid_seqs_lib.py
        # print(df_metadata.shape, df_metadata.columns)
        ndxs = filter_hs_by_date_new(df_metadata, _date_inf, _date_sup, date_mask="%d/%m/%Y")
        # print(">>> ndxs", _date_inf, _date_sup, df_metadata.shape)
        # print(">>> ndxs", len(ndxs))
        mseqProt.seq_records = [mseqProt.seq_records[i] for i in ndxs]
        mseqProt.seqs        = [mseqProt.seqs[i]        for i in ndxs]
        nSamples  = len(mseqProt.seqs)
        print(">>> seqs", nSamples)
        #--- must recalculate all entropy again !
        ret = shan.fast_calcHSannon_bias_correction_1pos(mseqProt, df_metadata, type="Protein", verbose=False)

        matPiDic = shan.dicPiList
        matNDic  = shan.dicNList
        hs       = shan.HShannonList
        maxVal   = np.max(hs)

        if increment == -1:
            end = -1
        else:
            end = ini+increment

        maximus = []

        maxi = len(hs)
        if end == -1:
            end = maxi
        elif end > maxi:
            end = maxi

        matPiDic = matPiDic[ini : end]
        matNDic  = matNDic[ini : end]
        hs       = hs[ini : end]
        maximus.append(np.max(hs))

        labelx = []; _seqx = []; _seqn = []; count = 1
        for jj in range(len(matPiDic)):
            #---- percent
            dic2 = matPiDic[jj]
            val = "; ".join(["%s %.2f%%"%(k, dic2[k]*100) for k in dic2.keys() if dic2[k] != 0])
            _seqx.append(val)

            valks = []
            for k in dic2.keys():
                if dic2[k] == 0: continue

                if dic2[k] == 1:
                    valks = [k]
                    break

                valks.append("%s %.2f%%"%(k, dic2[k]*100) )

            labelx.append(str(count) + '-' + "; ".join(valks))

            #---- n or #
            dic3 = matNDic[jj]
            val = "; ".join(["%s %d"%(k, dic3[k]) for k in dic3.keys() if dic3[k] != 0])
            _seqn.append(val)

            count += 1

        seqx = np.arange(ini, end)
        if i > 0:
            seqx = [x + delwidth*i for x in seqx]

        xmaxi = len(hs)
        maximus = np.max(maximus)
        dicdata = {'country' : [country]*len(hs), 'x': seqx, 'y': hs, 'aas': _seqx, 'nns': _seqn}

        return dicdata, nSamples, maxVal


    def on_button_submit(event):
        if txt_dummy.value == 'a':
            txt_dummy.value = 'b'
        else:
            txt_dummy.value = 'a'

    button_submit.on_event(ButtonClick, on_button_submit)

    #---------------- ojbects ----------------------------------------
    # min_month, max_month,

    controls = [group_dates, button_submit, group_protein, group_countries]
    # txt_date_inf, txt_date_sup, checkbox_countries
    for control in [txt_dummy, txt_alert]:
        if control == txt_alert:
            # control.on_change('value', lambda attr, old, new: update())
            pass
        else:
            control.on_change('value', lambda attr, old, new: update())

    inputs = column(*controls, width=int(plot_width/4), height=plot_height)
    inputs.sizing_mode = "fixed"
    lout = layout([ [deschtml], [inputs, p] ], sizing_mode="scale_both" )

    update()  # initial load of the data
    # curdoc().add_root(lout)

    return render_template('mutation_map.html', title=_('Mutation Map - SARS-CoV-2'), form=form, lout=lout)


@bp.route('/follow/<username>', methods=['POST'])
@login_required
def follow(username):
    form = EmptyForm()
    if form.validate_on_submit():
        user = User.query.filter_by(username=username).first()
        if user is None:
            flash(_('User %(username)s not found.', username=username))
            return redirect(url_for('main.index'))
        if user == current_user:
            flash(_('You cannot follow yourself!'))
            return redirect(url_for('main.user', username=username))
        current_user.follow(user)
        db.session.commit()
        flash(_('You are following %(username)s!', username=username))
        return redirect(url_for('main.user', username=username))
    else:
        return redirect(url_for('main.index'))


@bp.route('/unfollow/<username>', methods=['POST'])
@login_required
def unfollow(username):
    form = EmptyForm()
    if form.validate_on_submit():
        user = User.query.filter_by(username=username).first()
        if user is None:
            flash(_('User %(username)s not found.', username=username))
            return redirect(url_for('main.index'))
        if user == current_user:
            flash(_('You cannot unfollow yourself!'))
            return redirect(url_for('main.user', username=username))
        current_user.unfollow(user)
        db.session.commit()
        flash(_('You are not following %(username)s.', username=username))
        return redirect(url_for('main.user', username=username))
    else:
        return redirect(url_for('main.index'))


@bp.route('/translate', methods=['POST'])
@login_required
def translate_text():
    return jsonify({'text': translate(request.form['text'],
                                      request.form['source_language'],
                                      request.form['dest_language'])})

@bp.route('/search')
@login_required
def search():
    if not g.search_form.validate():
        return redirect(url_for('main.explore'))
    page = request.args.get('page', 1, type=int)
    posts, total = Post.search(g.search_form.q.data, page,
                               current_app.config['POSTS_PER_PAGE'])
    next_url = url_for('main.search', q=g.search_form.q.data, page=page + 1) \
        if total > page * current_app.config['POSTS_PER_PAGE'] else None
    prev_url = url_for('main.search', q=g.search_form.q.data, page=page - 1) \
        if page > 1 else None
    return render_template('search.html', title=_('Search'), posts=posts,
                           next_url=next_url, prev_url=prev_url)
