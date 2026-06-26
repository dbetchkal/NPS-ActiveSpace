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
        url = engine.url.render_as_string(hide_password=False)
        assert url == "postgresql://user%40domain:p%40ssword@10.0.0.1:5432/overflights"
