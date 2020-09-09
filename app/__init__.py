import logging
from logging.handlers import SMTPHandler, RotatingFileHandler
import os, sys, time, shutil

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from app.dashPlots.layout    import layout
from app.dashPlots.callbacks import register_callbacks

from flask import Flask, request, current_app
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from flask_login import LoginManager
from flask_mail import Mail
from flask_bootstrap import Bootstrap
from flask_moment import Moment
from flask_babel import Babel, lazy_gettext as _l
from flask.helpers import get_root_path
from flask_login import login_required
# from elasticsearch import Elasticsearch
from config import Config

sys.path.insert(1, '../src/')
from biobanco_lib   import *
from Sequence import *
from Shannon import *
from parallel import *
from covid_seqs_lib import *


db = SQLAlchemy()
migrate = Migrate()
login = LoginManager()
login.login_view = 'auth.login'
login.login_message = _l('Please log in to access this page.')
mail = Mail()
bootstrap = Bootstrap()
moment = Moment()
babel = Babel()
mydash = None

def create_app(config_class=Config):
    wantDash = True

    app = Flask(__name__)
    app.config.from_object(config_class)

    db.init_app(app)
    migrate.init_app(app, db)
    login.init_app(app)
    mail.init_app(app)
    bootstrap.init_app(app)
    moment.init_app(app)
    babel.init_app(app)

    app.file_index_css = rename_index_css()
    app.jinja_env.globals['file_index_css'] = app.file_index_css
    # print(">>> ", app.file_index_css)

    '''app.elasticsearch = Elasticsearch([app.config['ELASTICSEARCH_URL']]) \
        if app.config['ELASTICSEARCH_URL'] else None'''

    from app.errors import bp as errors_bp
    app.register_blueprint(errors_bp)

    from app.auth import bp as auth_bp
    app.register_blueprint(auth_bp, url_prefix='/auth')

    from app.main import bp as main_bp
    app.register_blueprint(main_bp)

    if not app.debug and not app.testing:
        if app.config['MAIL_SERVER']:
            auth = None
            if app.config['MAIL_USERNAME'] or app.config['MAIL_PASSWORD']:
                auth = (app.config['MAIL_USERNAME'],
                        app.config['MAIL_PASSWORD'])
            secure = None
            if app.config['MAIL_USE_TLS']:
                secure = ()
            mail_handler = SMTPHandler(
                mailhost=(app.config['MAIL_SERVER'], app.config['MAIL_PORT']),
                fromaddr='no-reply@' + app.config['MAIL_SERVER'],
                toaddrs=app.config['ADMINS'], subject='CENTD portal Failure',
                credentials=auth, secure=secure)
            mail_handler.setLevel(logging.ERROR)
            app.logger.addHandler(mail_handler)

        if app.config['LOG_TO_STDOUT']:
            stream_handler = logging.StreamHandler()
            stream_handler.setLevel(logging.INFO)
            app.logger.addHandler(stream_handler)
        else:
            if not os.path.exists('logs'):
                os.mkdir('logs')
            file_handler = RotatingFileHandler('logs/centd.log', maxBytes=10240, backupCount=10)
            file_handler.setFormatter(logging.Formatter(
                '%(asctime)s %(levelname)s: %(message)s '
                '[in %(pathname)s:%(lineno)d]'))
            file_handler.setLevel(logging.INFO)
            app.logger.addHandler(file_handler)

        app.logger.setLevel(logging.INFO)
        app.logger.info('CENTD portal startup')

    global mydash
    mydash = init_dashboard_first(app)
    # print("app __init__ mydash = ", mydash)

    return mydash

def rename_index_css():
    want_rename = True

    if want_rename:
        _root = 'app/static/css/'
        filename = "index.css"

        todels = [x for x in os.listdir(_root) if 'index_time_' in x]
        for fdel in todels:
            os.remove(os.path.join(_root, fdel))

        filenew = 'index_time_%s.css'%(time.ctime()).replace(' ','_').replace(":",'-')

        shutil.copy(os.path.join(_root, filename), os.path.join(_root,filenew) )
    else:
        filenew = 'index.css'

    return filenew


@babel.localeselector
def get_locale():
    return request.accept_languages.best_match(current_app.config['LANGUAGES'])


# https://medium.com/@olegkomarov_77860/how-to-embed-a-dash-app-into-an-existing-flask-app-ea05d7a2210b
def init_dashboard_first(app):

    # Meta tags for viewport responsiveness
    meta_viewport = {"name": "viewport", "content": "width=device-width, initial-scale=1, shrink-to-fit=no"}

    mydash = dash.Dash(__name__,
                         server=app,
                         url_base_pathname='/dashboard/',
                         assets_folder=get_root_path(__name__) + '/dashboard/assets/',
                         meta_tags=[meta_viewport])

    @mydash.callback(
        Output('heatmap-id', 'figure'),
        [Input('experiment-id', 'value')])
    def update_biobanc_results_plot(vid):
        # print(">>> update_graph", vid)
        return callback_biobanc_results_plot(vid)


    @mydash.callback(
        Output('entropic-plot-id', 'figure'),
        [Input('submit', 'n_clicks')],
        [State('country_picker', 'value'),
         State('protein_picker', 'value'),
         State('date_picker', 'start_date'),
         State('date_picker', 'end_date')])
    def update_entropic_plot(n_clicks, countryList, protein, start_date, end_date):
        return callback_update_entropic_plot(n_clicks, countryList, protein, start_date, end_date)

    with app.app_context():
        mydash.title = 'Dashapp 1'
        mydash.layout = layout
        register_callbacks(mydash)

    _protect_dashviews(mydash)

    return mydash

def callback_update_entropic_plot(n_clicks, countryList, protein, start_date, end_date):

    try:
        from covid_seqs_lib import mm
        print("####", mm)
    except:
        stri = "NOT FOUND - import mm"
        print(stri)
        return None

    if mm is None:
        stri = "NOT FOUND - mm"
        print(stri)
        return None

    # parallel (lib)
    matdic = run_multiprocess_countries(countryList, protein=protein,
                                        start_date=start_date, end_date=end_date,
                                        prjName=mm.prjName, isProtein=mm.isProtein,
                                        base_filename=mm.base_filename, root_data=mm.root_data,
                                        root_protein=mm.root_protein, root_figure=mm.root_protent_figure,
                                        root_entropy=mm.root_protent_entropy, root_html=mm.root_protent_html,
                                        root_templates=mm.root_templates, date_mask=mm.date_mask, cpus=mm.cpus)

    traces = []
    for imat in range(len(matdic)):
        country  = matdic[imat]['country']
        df       = matdic[imat]['df']
        nSamples = matdic[imat]['nSamples']
        maxVal   = matdic[imat]['maxVal']

        if len(df) == 0:
            continue

        traces.append(go.Bar(x = df.x, y= df.y, text=df.aas, name=country) )

    if traces == []:
        fig = {'data': [{'x': [1,2,3], 'y': [3,3,3,]}],
               'layout':{'title': "Nothing found"}}
    else:
        fig = {'data': traces,
               'layout':{
                     'title': {
                         'text': "Polymorphism/Mutations for '%s'"%(protein),
                         'font': {'family': "Arial",
                                  'size': 24,
                                  'color': "black" } },
                     'xaxis': { 'title': {
                         'text': "protein residues",
                         'font': {'family': "Arial",
                                  'size': 18,
                                  'color': "black" } } },
                     'yaxis': { 'title': {
                         'text': "Entropy - H(bits)",
                         'font': {'family': "Arial",
                                  'size': 18,
                                  'color': "black" } } },
                     'legend': {
                         'title': {
                             'text': "Countries",
                             'font': {'family': "Arial",
                                      'size': 18,
                                      'color': "black" } },
                        'font': {
                            'family': "Arial",
                            'size':   16,
                            'color': "black"}
                     },
                } }
    return fig


def callback_biobanc_results_plot(vid, symmetric=True, want_nSamples=False):
    # print(">>>> callback_biobanc_results_plot", vid)
    vid = int(vid)

    try:
        from biobanco_lib import phm
    except:
        stri = "NOT FOUND - phm"
        print(stri)
        return None

    eperiment = phm.experiment
    dicexp    = phm.dicexp
    width     = phm.width

    height          = dicexp[vid]['height']
    name            = dicexp[vid]['name']
    fontsize        = dicexp[vid]['fontsize']
    sufix           = dicexp[vid]['author']
    filename        = dicexp[vid]['filename']
    experiment_type = dicexp[vid]['experiment_type']
    the_control     = dicexp[vid]['control']

    print(">>>> callback_biobanc_results_plot eperiment:", eperiment, ' - ', name, "\n", filename)

    phm.bb = Biobanco(experiment_type, the_control, sufix, filename,
                      phm.root_data, phm.root_result, phm.exclude_exp,
                      nround=phm.nround, valpha=phm.valpha, verbose=phm.verbose)

    dfs, nSamples = phm.bb.create_heamap_table(verbose=False)

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

    template = 'plotly_white'; fontcolor='black'

    if want_nSamples:
        title = "LFC Heatmap: venoms x cytokines<br>%s model; cell type: '%s'<br>samples = %s"%(phm.bb.experiment_type, phm.bb.cell_type, nSamples)
    else:
        title = "LFC Heatmap: venoms x cytokines<br>%s model; cell type: '%s'"%(phm.bb.experiment_type, phm.bb.cell_type)

    return prepare_heatmap_fig(dfs, title, mini, maxi, width, height, template, fontsize, fontcolor)


def register_dashapps(app):
    from app.dashPlots.layout    import layout
    from app.dashPlots.callbacks import register_callbacks

    # Meta tags for viewport responsiveness
    meta_viewport = {"name": "viewport", "content": "width=device-width, initial-scale=1, shrink-to-fit=no"}

    global mydash
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


from app import models
