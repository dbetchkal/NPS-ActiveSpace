from nps_active_space.utils.helpers import create_overflights_engine


class TestCreateOverflightsEngine:
    def test_url_encodes_at_in_username_and_password(self):
        engine = create_overflights_engine(
            {
                "username": "user@domain",
                "password": "p@ssword",
                "host": "10.0.0.1",
                "port": "5432",
                "name": "overflights",
            }
        )
        assert engine.url.drivername == "postgresql+psycopg2"
        assert engine.url.username == "user@domain"
        assert engine.url.password == "p@ssword"
        assert engine.url.host == "10.0.0.1"
        assert engine.url.port == 5432
        assert engine.url.database == "overflights"
        url = engine.url.render_as_string(hide_password=False)
        assert "user%40domain" in url
        assert "p%40ssword" in url
