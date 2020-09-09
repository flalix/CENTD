import logging
from logging.handlers import SMTPHandler, RotatingFileHandler
import os, sys, time, shutil

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
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
    print("app __init__ mydash = ", mydash)
    mydash = init_dashboard_first(app)
    print("app __init__ mydash = ", mydash)

    '''
    global mydash
    mydash = init_dashboard_first(app)

    with app.app_context():
        # Import Dash application
        # from .plotlydash.dashboard import create_dashboard
        # dash_app = create_dashboard(app)
        global dboard
        dboard = Dashboard(app)
        return dboard.dash_app'''

    return mydash

def rename_index_css():
    want_rename = True

    if want_rename:
        _root = 'app/static/css/'
        filename = "index.css" # [x for x in os.listdir(_root) if 'index.css' in x][0]
        filenew = 'index_%s.css'%(time.ctime()).replace(' ','_').replace(":",'-')

        # print(">>> ", _root+filename, _root+filenew)
        shutil.copy(_root+filename, _root+filenew)
    else:
        filenew = 'index.css'

    return filenew


@babel.localeselector
def get_locale():
    return request.accept_languages.best_match(current_app.config['LANGUAGES'])


class Dashboard(object):
    def __init__(self, server):
        self.server   = server
        self.dash_app = self.create_dashboard(server)

    # https://hackersandslackers.com/plotly-dash-with-flask/
    def create_dashboard(self, server):
        """Create a Plotly Dash dashboard."""

        dash_app = dash.Dash(
            server=server,
            routes_pathname_prefix='/dashapp/',
            external_stylesheets=[
                '/static/css/styles.css',
            ]
        )

        # Create Dash Layout
        dash_app.layout = html.Div(id='dash-container')
        self.init_callbacks(dash_app)

        return dash_app

    def init_callbacks(self, dash_app):
        @dash_app.callback(Output('heatmap-id', 'figure'),
                          [Input('experiment-id', 'value')])
        def update_graph(vid):
            print("callback", vid)
            return vid

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
    def update_graph(vid):
        print(">>> update_graph", vid)
        return callback_update_graph(vid)

    with app.app_context():
        mydash.title = 'Dashapp 1'
        mydash.layout = layout
        register_callbacks(mydash)

    _protect_dashviews(mydash)

    return mydash

def callback_update_graph(vid, symmetric=True, want_nSamples=False):

    print(">>>> callback_update_graph", vid)
    vid = int(vid)

    try:
        from biobanco_lib import phm
    except:
        stri = "NOT FOUND - phm"
        print(stri)
        return stri

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

    print(">>>> callback_update_graph eperiment:", eperiment, ' - ', name, "\n", filename)

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
