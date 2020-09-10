import copy, os, random, sys, re, math
from os.path import dirname, join
import numpy as np
import pandas as pd
from collections import OrderedDict
from operator import itemgetter #, attrgetter
from datetime import datetime

from flask import render_template, flash, redirect, url_for, request, g, jsonify, current_app
from flask_login import current_user, login_required
from flask_babel import _, get_locale
from flask.helpers import get_root_path
from guess_language import guess_language

from app.dashPlots.layout    import layout
from app.dashPlots.callbacks import register_callbacks

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

from app import db
from app.main.forms import EditProfileForm, EmptyForm, PostForm , SearchForm, \
                           VirusForm, MutationMapForm, LabResultsForm
from app.models import User, Post
from app.translate import translate
from app.main import bp

from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import plotly.graph_objects as go
import plotly.express as pxydash.layout = layout
template = 'plotly_white'

#-- plot heatmap instance
phm = None

sys.path.insert(1, '../src/')
from Basic          import *
from Sequence       import *
from Shannon        import *
from covid_seqs_lib import *
from biobanco_lib   import *

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

experiments = ['Chondrocytes (Paola)', "Synoviocites (Isadora)", "Osteoclast (Eduardo)",
               'Inflammasome (Nancy and Olga)', 'Macrophages (Douglas and Fernanda)',
               'Neurons (Michelle)']

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
    return render_template('lab_results_202008.html', title='Screening results', experiments=experiments)

@bp.route('/show_exp_results/<string:experiment>')
@login_required
def show_exp_results(experiment):

    from app import mydash
    print(">>> experiment", experiment)
    # layout = prepare_scatter_plot()

    # from copy import deepcopy
    try:
        global phm
        phm, layout = define_biobanc_experiment(experiment)
    except:
        print(">>> NOT got layout ", experiment)
        pass
    mydash.layout = layout

    # print(url_for('/dashboard/'))
    return redirect(url_for('/dashboard/') )

    # return render_template('heatmap.html', title='Results for '+experiment, experiment=experiment)

'''@mydash.callback(Output('heatmap-id', 'figure'),
                [Input('experiment-id', 'value')])'''
def callback_update_graph(vid, symmetric=True, want_nSamples=False):

    print(">>>> callback", vid)
    try:
        vid = int(vid)
    except:
        vid = 1

    eperiment = phm.experiment
    dicexp    = phm.dicexp
    width     = phm.width

    height          = dicexp[vid]['height'];
    sufix           = dicexp[vid]['author']
    filename        = dicexp[vid]['filename']
    experiment_type = dicexp[vid]['experiment_type']
    the_control     = dicexp[vid]['control']

    phm.bb = Biobanco(phm.experiment_type, phm.the_control, phm.sufix, filename,
                      phm.root_data, phm.root_result, phm.exclude_exp,
                      nround=phm.nround, valpha=phm.valpha, verbose=phm.verbose)

    dfs, nSamples = phm.bb.create_heamap_table(verbose=verbose)

    if dfs is None:
        print("Nothing found")
        return None

    maxi = dfs.lfc.max()
    mini = dfs.lfc.min()

    if not symmetric:
        if maxi < _maxi: maxi = _maxi
        if mini > _mini: mini = _mini
    else:
        if abs(mini) > maxi:
            maxi = abs(mini)
        else:
            mini = -maxi

    if want_nSamples:
        title = "LFC Heatmap: venoms x cytokines<br>%s model; cell type: '%s'<br>samples = %s"%(phm.bb.experiment_type, phm.bb.cell_type, nSamples)
    else:
        title = "LFC Heatmap: venoms x cytokines<br>%s model; cell type: '%s'"%(phm.bb.experiment_type, phm.bb.cell_type)

    return prepare_heatmap_fig(dfs, title, mini, maxi, width, height, template, fontsize, fontcolor)


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

    print(">>> mutation_map")
    # layout = prepare_scatter_plot()
    layout = init_mutation_layout(cpus=6)

    # print(url_for('/dashboard/'))
    return redirect(url_for('/dashboard/') )


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
