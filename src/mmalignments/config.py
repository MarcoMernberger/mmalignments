from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    environment: str = "development"
    db_url: str | None = None
    redis_url: str | None = None

    class Config:
        env_file = ".env"


settings = Settings()
